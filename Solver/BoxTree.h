#ifndef BOXTREE_H
#define BOXTREE_H

#include <vector>
#include <limits>
#include <map>
#include <set>
#include <iostream>
#include <numeric>

#include <Eigen/Dense>

#include <omp.h>
#include "mpi.h"
#include "mkl.h"

#include "Level.h"
#include "Vector.h"
#include "parameters.h"
#include "Interpolator.h"

#include <sys/time.h>
#include <sys/resource.h>

#include <fstream>

class Level;
class BoxTree;

const Interpolator<Chebyshev> IPSCHEME_;

class BoxTree 
{

    private:

        int comm_rank_;
        int comm_size_;
        MPI_Comm mpi_comm_;

        std::vector<long long> split_points_orig_;
        std::vector<int> recv_counts_orig_;
        std::vector<int> displs_orig_;
        long long N_loc_orig_;
        long long point_low_;
        long long point_up_;

        std::vector<double> x_;
        std::vector<double> y_;
        std::vector<double> z_;

        std::vector<double> normal_x_;
        std::vector<double> normal_y_;
        std::vector<double> normal_z_;        

        std::vector<long long> sorting_;

        int nlevels_;
        std::complex<double> wavenumber_;  

        std::complex<double> coupling_parameter_;
        std::unordered_map<long long, std::unordered_set<long long>> precomputations_data_;

        std::vector<Level*> levels_;

        std::vector<std::complex<double>> solution_;
        
        int P_;

        std::vector<int> n_pts_per_patch_;

    private:

        void Initialize() 
        {

            InitializeMPI();

            InitializeBoxesAndLevels();

            if (USE_HIGH_ORDER) {

                InitializeRelevantConeSegmentsHO(); 

            } else {

                InitializeRelevantConeSegments();

            }

            InitializeSolutionAndConeSegmentCoefficients();

        }

        void InitializeMPI()
        {

            int init;
            MPI_Initialized(&init);

            if (!init)
                throw std::logic_error("MPI needs to be initialized before calling the BoxTree\n");

            MPI_Comm_size(mpi_comm_, &comm_size_);
            MPI_Comm_rank(mpi_comm_, &comm_rank_);

            point_low_ = split_points_orig_[comm_rank_];
            point_up_ = split_points_orig_[comm_rank_+1];

        }

        void ComputeBB(std::array<double, 3>& min, std::array<double ,3>& max) const
        {

            std::array<double, 3> min_loc = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
            std::array<double, 3> max_loc = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};

            for (long long i = 0; i < N_loc_orig_; i++) {

                if (x_[i] < min_loc[0]) min_loc[0] = x_[i];
                if (y_[i] < min_loc[1]) min_loc[1] = y_[i];
                if (z_[i] < min_loc[2]) min_loc[2] = z_[i];
                
                if (x_[i] > max_loc[0]) max_loc[0] = x_[i];
                if (y_[i] > max_loc[1]) max_loc[1] = y_[i];
                if (z_[i] > max_loc[2]) max_loc[2] = z_[i];
                
            }

            MPI_Allreduce(&min_loc[0], &min[0], 3, MPI_DOUBLE, MPI_MIN, mpi_comm_);
            MPI_Allreduce(&max_loc[0], &max[0], 3, MPI_DOUBLE, MPI_MAX, mpi_comm_);

            max[0] += 1e-10;
            max[1] += 1e-10;
            max[2] += 1e-10;

            min[0] -= 1e-10;
            min[1] -= 1e-10;
            min[2] -= 1e-10;
            
        }

        void ComputeNLevels(const std::array<double, 3>& min, const double boxsize)
        {

            int nlevels = 0;

            bool split = true;

            while (split) {

                nlevels++;

                const long long nboxes_axis = 1LL << (nlevels - 1);
                
                std::unordered_map<long long, long long> morton_code_counts_loc;                

                for (long long i = 0; i < N_loc_orig_; i++) {

                    const long long morton_box = Level::Point2Morton(x_[i], y_[i], z_[i], min[0], min[1], min[2], boxsize / nboxes_axis, nlevels - 1);
 
                    morton_code_counts_loc[morton_box]++;                  

                }

                std::vector<long long> morton_box_loc;
                std::vector<long long> morton_box_count_loc;

                morton_box_loc.reserve(morton_code_counts_loc.size());
                morton_box_count_loc.reserve(morton_code_counts_loc.size());

                for (const auto & [morton_box, count] : morton_code_counts_loc) {

                    morton_box_loc.push_back(morton_box);
                    morton_box_count_loc.push_back(count);

                }

                morton_code_counts_loc.clear();

                MPI_Count size_loc = morton_box_loc.size();

                std::vector<MPI_Count> recv_counts(comm_size_);
                std::vector<MPI_Aint> displs(comm_size_);

                MPI_Allgather_c(&size_loc, 1, MPI_COUNT, &recv_counts[0], 1, MPI_COUNT, mpi_comm_);

                MPI_Count size_all = 0;                

                for (int rank = 0; rank < comm_size_; rank++) {

                    displs[rank] = size_all;
                    size_all += recv_counts[rank];

                }

                std::vector<long long> morton_box_all(static_cast<size_t>(size_all));
                std::vector<long long> morton_box_count_all(static_cast<size_t>(size_all));

                MPI_Allgatherv_c(&morton_box_loc[0], size_loc, MPI_LONG_LONG, &morton_box_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
                MPI_Allgatherv_c(&morton_box_count_loc[0], size_loc, MPI_LONG_LONG, &morton_box_count_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);

                std::vector<long long>().swap(morton_box_loc);
                std::vector<long long>().swap(morton_box_count_loc);

                std::unordered_map<long long, long long> morton_code_counts_all;

                for (MPI_Count i = 0; i < size_all; i++) {

                    const long long morton_box = morton_box_all[static_cast<size_t>(i)];
                    const long long counts = morton_box_count_all[static_cast<size_t>(i)];

                    morton_code_counts_all[morton_box] += counts;

                }

                std::vector<long long>().swap(morton_box_all);
                std::vector<long long>().swap(morton_box_count_all);

                long long max_elem = 0;

                for (const auto & [morton_box, count] : morton_code_counts_all) {

                    if (max_elem < count) max_elem = count;

                }

                split = max_elem > MAX_ELEMS_LEAF;

                morton_code_counts_all.clear();

            }

            nlevels_ = nlevels;

        }

        void SortBox(const std::array<double, 3>& min, const double boxsize) 
        {

            std::vector<long long> morton_code(N_loc_orig_);

            for (long long i = 0; i < N_loc_orig_; i++) {

                const long long morton = Level::Point2Morton(x_[i], y_[i], z_[i], min[0], min[1], min[2], boxsize, nlevels_-1);
                morton_code[i] = morton;

            }                       

            std::vector<long long> sorting_all;
            std::vector<long long> morton_code_all;

            std::vector<double> xyz_all;
            std::vector<double> normal_xyz_all;

            std::vector<double> tmp_xyz_all;
            std::vector<double> tmp_normal_xyz_all;

            long long N = split_points_orig_[comm_size_];

            if (comm_rank_ == 0) {

                sorting_all.resize(N);
                morton_code_all.resize(N);
                
                xyz_all.resize(N);
                normal_xyz_all.resize(N);

                tmp_xyz_all.resize(N);
                tmp_normal_xyz_all.resize(N);

            }

            MPI_Gatherv(&sorting_[0], N_loc_orig_, MPI_LONG_LONG, &sorting_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_LONG_LONG, 0, mpi_comm_);
            MPI_Gatherv(&morton_code[0], N_loc_orig_, MPI_LONG_LONG, &morton_code_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_LONG_LONG, 0, mpi_comm_);

            if (comm_rank_ == 0) {

                new_sort(sorting_all, morton_code_all);

            }               

            std::vector<long long>().swap(morton_code);
            std::vector<long long>().swap(morton_code_all);

            MPI_Gatherv(&x_[0], N_loc_orig_, MPI_DOUBLE, &xyz_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, 0, mpi_comm_);
            MPI_Gatherv(&normal_x_[0], N_loc_orig_, MPI_DOUBLE, &normal_xyz_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, 0, mpi_comm_);
 
            if (comm_rank_ == 0) {                

                for (long long i = 0; i < N; i++) {

                    tmp_xyz_all[i] = xyz_all[sorting_all[i]];
                    tmp_normal_xyz_all[i] = normal_xyz_all[sorting_all[i]];

                }

            }

            MPI_Scatterv(&tmp_xyz_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, &x_[0], N_loc_orig_, MPI_DOUBLE, 0, mpi_comm_);
            MPI_Scatterv(&tmp_normal_xyz_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, &normal_x_[0], N_loc_orig_, MPI_DOUBLE, 0, mpi_comm_);
            
            MPI_Gatherv(&y_[0], N_loc_orig_, MPI_DOUBLE, &xyz_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, 0, mpi_comm_);
            MPI_Gatherv(&normal_y_[0], N_loc_orig_, MPI_DOUBLE, &normal_xyz_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, 0, mpi_comm_);
 
            if (comm_rank_ == 0) {                

                for (long long i = 0; i < N; i++) {

                    tmp_xyz_all[i] = xyz_all[sorting_all[i]];
                    tmp_normal_xyz_all[i] = normal_xyz_all[sorting_all[i]];

                }

            }

            MPI_Scatterv(&tmp_xyz_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, &y_[0], N_loc_orig_, MPI_DOUBLE, 0, mpi_comm_);
            MPI_Scatterv(&tmp_normal_xyz_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, &normal_y_[0], N_loc_orig_, MPI_DOUBLE, 0, mpi_comm_);

            MPI_Gatherv(&z_[0], N_loc_orig_, MPI_DOUBLE, &xyz_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, 0, mpi_comm_);
            MPI_Gatherv(&normal_z_[0], N_loc_orig_, MPI_DOUBLE, &normal_xyz_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, 0, mpi_comm_);
 
            if (comm_rank_ == 0) {                

                for (long long i = 0; i < N; i++) {

                    tmp_xyz_all[i] = xyz_all[sorting_all[i]];
                    tmp_normal_xyz_all[i] = normal_xyz_all[sorting_all[i]];

                }

            }

            MPI_Scatterv(&tmp_xyz_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, &z_[0], N_loc_orig_, MPI_DOUBLE, 0, mpi_comm_);
            MPI_Scatterv(&tmp_normal_xyz_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, &normal_z_[0], N_loc_orig_, MPI_DOUBLE, 0, mpi_comm_);
            
            std::vector<double>().swap(xyz_all);
            std::vector<double>().swap(normal_xyz_all);
            std::vector<double>().swap(tmp_xyz_all);
            std::vector<double>().swap(tmp_normal_xyz_all);

            MPI_Scatterv(&sorting_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_LONG_LONG, &sorting_[0], N_loc_orig_, MPI_LONG_LONG, 0, mpi_comm_);

            std::vector<long long>().swap(sorting_all);

            std::vector<double> x_all(N);
            std::vector<double> y_all(N);
            std::vector<double> z_all(N);

            std::vector<double> normal_x_all(N);
            std::vector<double> normal_y_all(N);
            std::vector<double> normal_z_all(N);

            MPI_Allgatherv(&x_[0], N_loc_orig_, MPI_DOUBLE, &x_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, mpi_comm_);
            MPI_Allgatherv(&y_[0], N_loc_orig_, MPI_DOUBLE, &y_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, mpi_comm_);
            MPI_Allgatherv(&z_[0], N_loc_orig_, MPI_DOUBLE, &z_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, mpi_comm_);

            MPI_Allgatherv(&normal_x_[0], N_loc_orig_, MPI_DOUBLE, &normal_x_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, mpi_comm_);
            MPI_Allgatherv(&normal_y_[0], N_loc_orig_, MPI_DOUBLE, &normal_y_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, mpi_comm_);
            MPI_Allgatherv(&normal_z_[0], N_loc_orig_, MPI_DOUBLE, &normal_z_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE, mpi_comm_);

            x_ = x_all;
            std::vector<double>().swap(x_all);
            y_ = y_all;
            std::vector<double>().swap(y_all);
            z_ = z_all;
            std::vector<double>().swap(z_all);

            normal_x_ = normal_x_all;
            std::vector<double>().swap(normal_x_all);
            normal_y_ = normal_y_all;
            std::vector<double>().swap(normal_y_all);
            normal_z_ = normal_z_all;
            std::vector<double>().swap(normal_z_all);

            sorting_all.resize(N);

            MPI_Allgatherv(&sorting_[0], N_loc_orig_, MPI_LONG_LONG, &sorting_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_LONG_LONG, mpi_comm_);

            sorting_ = sorting_all;
            std::vector<long long>().swap(sorting_all);

        }

        void InitializeLevelDBoxesData() 
        {

            long long old_morton = -1;
            long long npoints;

            for (long long i = 0; i < N_loc_orig_; i++) {

                const long long idx = i + point_low_;

                const long long morton = levels_.back()->Point2Morton(x_[idx], y_[idx], z_[idx]);

                if (morton != old_morton) {

                    if (i != 0) {

                        levels_.back()->mortonbox2discretizationpoints_[old_morton] = {idx - npoints, npoints};

                    }

                    old_morton = morton;
                    npoints = 1;
                    levels_.back()->mortonidofrelboxes_.insert(morton);

                } else {

                    npoints++;

                }

            }

            levels_.back()->mortonbox2discretizationpoints_[old_morton] = {point_up_ - npoints, npoints}; 

        }

        void InitializeRelBoxesAllLevels() 
        {

            std::unordered_set<long long> relevantparentmortonown;
            std::unordered_set<long long> relevantparentmortonothers;
            std::unordered_set<long long> relevantmortonown;
            std::unordered_set<long long> relevantmortonothers;

            std::vector<long long> relevantmorton_loc(levels_.back()->mortonidofrelboxes_.begin(), levels_.back()->mortonidofrelboxes_.end());
            int numrelboxes_loc = relevantmorton_loc.size();

            std::vector<int> numrelboxes_all(comm_size_);
            MPI_Allgather(&numrelboxes_loc, 1, MPI_INT, &numrelboxes_all[0], 1, MPI_INT, mpi_comm_);

            std::vector<int> displs(comm_size_, 0);
            long long numrelboxes_total = 0;

            for (int i = 0; i < comm_size_; i++) {

                numrelboxes_total += numrelboxes_all[i];

                if (i != 0)
                    displs[i] = displs[i-1] + numrelboxes_all[i-1];

            }
            
            std::vector<long long> mortonidofrelboxes_all(numrelboxes_total);
            MPI_Allgatherv(&relevantmorton_loc[0], numrelboxes_loc, MPI_LONG_LONG, &mortonidofrelboxes_all[0], &numrelboxes_all[0], &displs[0], MPI_LONG_LONG, mpi_comm_);

            for (long long i = 0; i < mortonidofrelboxes_all.size(); i++) {

                if (levels_.back()->mortonidofrelboxes_.count(mortonidofrelboxes_all[i]) != 0) {
                    relevantmortonown.insert(mortonidofrelboxes_all[i]);
                } else {
                    relevantmortonothers.insert(mortonidofrelboxes_all[i]);
                }

            }

            std::vector<long long>().swap(relevantmorton_loc);
            std::vector<long long>().swap(mortonidofrelboxes_all);

            for (int i = nlevels_-1; i >= 2; i--) {

                if (i != nlevels_-1) {

                    levels_[i]->mortonidofrelboxes_.insert(relevantmortonown.begin(), relevantmortonown.end());

                }

                for (const auto & morton_box : relevantmortonown) {

                    std::vector<long long> neighbors_and_cousins = levels_[i]->GetAllNeighboursAndCousins(morton_box);

                    for (long long k = 0; k < neighbors_and_cousins.size(); k++) {

                        const long long nc_morton_box = neighbors_and_cousins[k];

                        if (relevantmortonown.count(nc_morton_box) != 0 || relevantmortonothers.count(nc_morton_box) != 0) {
                            levels_[i]->mortonidofrelboxes_2_.insert(nc_morton_box);
                        }

                    }

                    relevantparentmortonown.insert(static_cast<long long>(morton_box/8));

                }

                for (const auto & morton_box : relevantmortonothers) {

                    const long long parent_morton = static_cast<long long>(morton_box/8);

                    if (relevantparentmortonown.count(parent_morton) == 0) {

                        relevantparentmortonothers.insert(parent_morton);

                    }

                }

                std::swap(relevantmortonown, relevantparentmortonown);
                relevantparentmortonown.clear();
                std::swap(relevantmortonothers, relevantparentmortonothers);
                relevantparentmortonothers.clear();

            }         
            
        }

        void InitializeBoxesDataAllLevels() 
        {           

            std::vector<long long> morton_box_loc;
            std::vector<long long> position_loc;
            std::vector<long long> size_loc;

            for (const auto & [morton_box, values] : levels_.back()->mortonbox2discretizationpoints_) {

                morton_box_loc.push_back(morton_box);
                position_loc.push_back(values[0]);
                size_loc.push_back(values[1]);

            }

            int total_size_loc = morton_box_loc.size();

            std::vector<int> recv_counts(comm_size_);
            MPI_Allgather(&total_size_loc, 1, MPI_INT, &recv_counts[0], 1, MPI_INT, mpi_comm_);

            std::vector<int> displs(comm_size_, 0);
            long long total_size = 0;

            for (int i = 0; i < comm_size_; i++) {

                total_size += recv_counts[i];

                if (i != 0) {
                    displs[i] = displs[i-1] + recv_counts[i-1];
                }

            }

            std::vector<long long> morton_box_all(total_size);
            std::vector<long long> position_all(total_size);
            std::vector<long long> size_all(total_size);

            MPI_Allgatherv(&morton_box_loc[0], total_size_loc, MPI_LONG_LONG, &morton_box_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
            MPI_Allgatherv(&position_loc[0], total_size_loc, MPI_LONG_LONG, &position_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
            MPI_Allgatherv(&size_loc[0], total_size_loc, MPI_LONG_LONG, &size_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);

            std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_all;

            for (long long i = 0; i < total_size; i++) {

                const long long morton_box = morton_box_all[i];
                const long long position = position_all[i];
                const long long size = size_all[i];

                if (mortonbox2discretizationpoints_all.count(morton_box) == 0) {

                    mortonbox2discretizationpoints_all[morton_box] = {position, size};

                } else {

                    const std::array<long long, 2> actual_value = mortonbox2discretizationpoints_all[morton_box];

                    mortonbox2discretizationpoints_all[morton_box] = {std::min(actual_value[0], position), actual_value[1] + size};

                }

            }            

            std::vector<long long>().swap(morton_box_loc);
            std::vector<long long>().swap(position_loc);
            std::vector<long long>().swap(size_loc);

            std::vector<long long>().swap(morton_box_all);
            std::vector<long long>().swap(position_all);
            std::vector<long long>().swap(size_all); 

            std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_loc = levels_.back()->mortonbox2discretizationpoints_;
            
            for (int i = nlevels_-1; i >= 2; i--) {

                if (i != nlevels_-1) {

                    levels_[i]->mortonbox2discretizationpoints_ = mortonbox2discretizationpoints_loc;

                }

                for (const auto & morton_box : levels_[i]->mortonidofrelboxes_2_) {

                    levels_[i]->mortonbox2discretizationpoints_2_[morton_box] = mortonbox2discretizationpoints_all[morton_box];

                }

                std::unordered_map<long long, std::array<long long, 2>> new_mortonbox2discretizationpoints_all;
                std::unordered_map<long long, std::array<long long, 2>> new_mortonbox2discretizationpoints_loc;

                for (const auto [morton_box, values] : mortonbox2discretizationpoints_all) {

                    const long long new_morton_box = static_cast<long long>(morton_box/8);

                    if (new_mortonbox2discretizationpoints_all.count(new_morton_box) == 0) {

                        new_mortonbox2discretizationpoints_all[new_morton_box] = values;

                    } else {

                        new_mortonbox2discretizationpoints_all[new_morton_box][0] = std::min(new_mortonbox2discretizationpoints_all[new_morton_box][0], values[0]);
                        new_mortonbox2discretizationpoints_all[new_morton_box][1] += values[1];

                    }

                }

                mortonbox2discretizationpoints_all = new_mortonbox2discretizationpoints_all;
                std::unordered_map<long long, std::array<long long, 2>>{}.swap(new_mortonbox2discretizationpoints_all);

                for (const auto [morton_box, values] : mortonbox2discretizationpoints_loc) {

                    const long long new_morton_box = static_cast<long long>(morton_box/8);

                    if (new_mortonbox2discretizationpoints_loc.count(new_morton_box) == 0) {

                        new_mortonbox2discretizationpoints_loc[new_morton_box] = values;

                    } else {

                        new_mortonbox2discretizationpoints_loc[new_morton_box][0] = std::min(new_mortonbox2discretizationpoints_loc[new_morton_box][0], values[0]);
                        new_mortonbox2discretizationpoints_loc[new_morton_box][1] += values[1];

                    }

                }

                mortonbox2discretizationpoints_loc = new_mortonbox2discretizationpoints_loc;
                std::unordered_map<long long, std::array<long long, 2>>{}.swap(new_mortonbox2discretizationpoints_loc);

            }        

        }

        void InitializeBoxesAndLevels() 
        {            

            std::array<double, 3> min, max;

            ComputeBB(min, max);
            
            const double boxsize = std::max(max[0] - min[0], std::max(max[1] - min[1], max[2] - min[2]));

            sorting_ = std::vector<long long>(N_loc_orig_);

            for (long long i = 0; i < N_loc_orig_; i++) {

                sorting_[i] = i + point_low_;

            }

            if ((nlevels_ == -1) && USE_ADAPTIVITY) {

                ComputeNLevels(min, boxsize);

            }

            SortBox(min, boxsize / (1 << (nlevels_ - 1)));
            MPI_Barrier(mpi_comm_);

            levels_.resize(nlevels_, nullptr);

            #pragma omp parallel for
            for (int i = 0; i < nlevels_; i++) {

                levels_[i] = new Level(i, min[0], min[1], min[2], boxsize/(1 << i), wavenumber_);

            }
            MPI_Barrier(mpi_comm_);

            InitializeLevelDBoxesData();
            MPI_Barrier(mpi_comm_);

            InitializeRelBoxesAllLevels();
            MPI_Barrier(mpi_comm_);

            InitializeBoxesDataAllLevels();
            MPI_Barrier(mpi_comm_);    

        }

        void GetRelevantConeSegmentsDueToCousinSurface(const int level, 
            std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment)
        {

            std::unordered_map<long long, std::unordered_set<long long>> tmpmortonboxnonrelconesegmentall;

            #pragma omp parallel
            {

            long long old_morton_box = -1;
            std::vector<long long> cousins;
            std::unordered_map<long long, std::unordered_set<long long>> tmpmortonboxnonrelconesegment; 

            #pragma omp for
            for (long long i = point_low_; i < point_up_; i++) {

                const double x = x_[i];
                const double y = y_[i];
                const double z = z_[i];

                const long long morton_box = levels_[level]->Point2Morton(x, y, z);

                if (morton_box != old_morton_box) {

                    cousins = levels_[level]->GetCousins(morton_box);
                    old_morton_box = morton_box;

                }

                for (int j = 0; j < cousins.size(); j++) {

                    const long long morton_cousin = cousins[j];

                    const long long size_cousin = levels_[level]->mortonbox2discretizationpoints_2_[morton_cousin][1];

                    if (USE_ADAPTIVITY && (size_cousin <= MAX_ELEMS_LEAF)) {

                        continue;

                    }

                    const long long nonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_cousin, x, y, z);

                    tmpmortonboxnonrelconesegment[morton_cousin].insert(nonrel_conesegment);

                }

            }

            for (auto iter = tmpmortonboxnonrelconesegment.begin(); iter != tmpmortonboxnonrelconesegment.end(); iter++) {

                #pragma omp critical
                tmpmortonboxnonrelconesegmentall[iter->first].insert(iter->second.begin(), iter->second.end());

            }
            
            }

            for (auto i = tmpmortonboxnonrelconesegmentall.begin(); i != tmpmortonboxnonrelconesegmentall.end(); i++) {

                std::vector<long long> tmp(i->second.begin(), i->second.end());
                std::sort(tmp.begin(), tmp.end());
                mortonboxnonrelconesegment[i->first] = tmp;

            }

        }  

        void GetRelevantConeSegmentsDueToCousinSurface(const int level, 
            std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment,
            std::unordered_map<long long, std::unordered_map<long long, double>> & min_radius_conesegment,
            std::unordered_map<long long, std::unordered_map<long long, double>> & max_radius_conesegment)
        {

            std::unordered_map<long long, std::unordered_set<long long>> tmpmortonboxnonrelconesegmentall;

            #pragma omp parallel
            {

            long long old_morton_box = -1;
            std::vector<long long> cousins;
            std::unordered_map<long long, std::unordered_set<long long>> tmpmortonboxnonrelconesegment; 
            std::unordered_map<long long, std::unordered_map<long long, double>> tmp_min_radius_conesegment;
            std::unordered_map<long long, std::unordered_map<long long, double>> tmp_max_radius_conesegment;            

            #pragma omp for
            for (long long i = point_low_; i < point_up_; i++) {

                const double x = x_[i];
                const double y = y_[i];
                const double z = z_[i];

                const long long morton_box = levels_[level]->Point2Morton(x, y, z);

                if (morton_box != old_morton_box) {

                    cousins = levels_[level]->GetCousins(morton_box);
                    old_morton_box = morton_box;

                }

                for (int j = 0; j < cousins.size(); j++) {

                    const long long morton_cousin = cousins[j];

                    const long long size_cousin = levels_[level]->mortonbox2discretizationpoints_2_[morton_cousin][1];

                    if (USE_ADAPTIVITY && (size_cousin <= MAX_ELEMS_LEAF)) {

                        continue;

                    }

                    const long long nonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_cousin, x, y, z);

                    tmpmortonboxnonrelconesegment[morton_cousin].insert(nonrel_conesegment);

                    std::array<double, 3> center;
                    levels_[level]->GetBoxCenter(morton_cousin, center[0], center[1], center[2]);

                    const double diff = std::sqrt((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2]));

                    if ((tmp_min_radius_conesegment.count(morton_cousin) == 0) || (tmp_min_radius_conesegment.at(morton_cousin).count(nonrel_conesegment) == 0)) {

                        tmp_min_radius_conesegment[morton_cousin][nonrel_conesegment] = diff - 1E-10;

                    } else {

                        tmp_min_radius_conesegment[morton_cousin][nonrel_conesegment] = std::min({tmp_min_radius_conesegment[morton_cousin][nonrel_conesegment], diff - 1E-10});

                    }

                    if ((tmp_max_radius_conesegment.count(morton_cousin) == 0) || (tmp_max_radius_conesegment.at(morton_cousin).count(nonrel_conesegment) == 0)) {

                        tmp_max_radius_conesegment[morton_cousin][nonrel_conesegment] = diff + 1E-10;

                    } else {

                        tmp_max_radius_conesegment[morton_cousin][nonrel_conesegment] = std::max({tmp_max_radius_conesegment[morton_cousin][nonrel_conesegment], diff + 1E-10});

                    }

                }

            }

            for (auto iter = tmpmortonboxnonrelconesegment.begin(); iter != tmpmortonboxnonrelconesegment.end(); iter++) {

                #pragma omp critical
                tmpmortonboxnonrelconesegmentall[iter->first].insert(iter->second.begin(), iter->second.end());

            }

            for (const auto & [key_1, map_1] : tmp_min_radius_conesegment) {
                for (const auto & [key_2, value] : map_1) {

                    #pragma omp critical
                    if ((min_radius_conesegment.count(key_1) == 0) || (min_radius_conesegment.at(key_1).count(key_2) == 0)) {

                        min_radius_conesegment[key_1][key_2] = value;

                    } else {

                        min_radius_conesegment[key_1][key_2] = std::min({min_radius_conesegment[key_1][key_2], value});

                    }

                }
            }

            for (const auto & [key_1, map_1] : tmp_max_radius_conesegment) {
                for (const auto & [key_2, value] : map_1) {

                    #pragma omp critical
                    if ((max_radius_conesegment.count(key_1) == 0) || (max_radius_conesegment.at(key_1).count(key_2) == 0)) {

                        max_radius_conesegment[key_1][key_2] = value;

                    } else {

                        max_radius_conesegment[key_1][key_2] = std::max({max_radius_conesegment[key_1][key_2], value});

                    }

                }
            }
            
            }

            for (auto i = tmpmortonboxnonrelconesegmentall.begin(); i != tmpmortonboxnonrelconesegmentall.end(); i++) {

                std::vector<long long> tmp(i->second.begin(), i->second.end());
                std::sort(tmp.begin(), tmp.end());
                mortonboxnonrelconesegment[i->first] = tmp;

            }

        }
        
        std::vector<long long> GetMortonRelChildren(const int level, const long long morton_box) const 
        {

            std::vector<long long> children;

            if (level >= nlevels_ - 1)
                return children;

            children.reserve(8);

            for (long long morton_children = morton_box*8; morton_children < (morton_box+1)*8; morton_children++) {

                if (levels_[level+1]->IsRelevant(morton_children))
                    children.push_back(morton_children);

            }

            return children;

        }        

        void GetRelevantConeSegmentsDueToParents(const int level, 
            std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment) 
        {

            std::unordered_map<long long, std::unordered_set<long long>> mortonboxnonrelconesegment_all;

            #pragma omp parallel
            {

            double x, y, z;
            long long old_cocentered_morton_box = -1;
            std::vector<long long> morton_children;
            std::unordered_map<long long, std::array<std::unordered_map<long long, std::vector<int>>, 8>> protobox;
            std::unordered_map<long long, std::unordered_set<long long>> mortonboxnonrelconesegment_thread;

            #pragma omp for
            for (long long iter = 0; iter < levels_[level-1]->split_points_relconesegments_[comm_rank_+1]-levels_[level-1]->split_points_relconesegments_[comm_rank_]; iter++) {

                const long long cocentered_morton_box = levels_[level-1]->relconesegment2cocenteredmortonboxid_[iter];
                const long long nonrel_conesegment = levels_[level-1]->relconesegment2nonrelconesegment_[iter];

                if (cocentered_morton_box != old_cocentered_morton_box) {

                    morton_children = GetMortonRelChildren(level-1, cocentered_morton_box);

                    if (morton_children.size() == 0) {
                        throw std::logic_error("Every cocentered morton box must have children in get relevant cone segments function.");
                    }

                    old_cocentered_morton_box = cocentered_morton_box;

                }

                for (int childiter = 0; childiter < morton_children.size(); childiter++) {

                    const long long morton_child = morton_children[childiter];
                    const long long child_pos = morton_child % 8;

                    if (protobox.count(nonrel_conesegment) == 0 || protobox.at(nonrel_conesegment)[child_pos].size() == 0) {

                        std::unordered_map<long long, std::vector<int>> & tmpchildproto = protobox[nonrel_conesegment][child_pos];

                        for (int interpiter = 0; interpiter < IPSCHEME_.GetNInterpPoints(); interpiter++) {

                            IPSCHEME_.GetInterpolationPoint(interpiter, x, y, z);

                            levels_[level-1]->Cheb2Cart(cocentered_morton_box, nonrel_conesegment, x, y, z);

                            long long locnonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_child, x, y, z);

                            tmpchildproto[locnonrel_conesegment].push_back(interpiter);

                        }

                    }

                }

            }

            #pragma omp critical
            {

            for (auto iter = protobox.begin(); iter != protobox.end(); iter++) {

                if (levels_[level]->protobox_.count(iter->first) == 0) {

                    levels_[level]->protobox_[iter->first] = iter->second;

                } else {

                    for (int childiter = 0; childiter < 8; childiter++) {

                        if (iter->second[childiter].size() != 0 && levels_[level]->protobox_[iter->first][childiter].size() == 0) {

                            levels_[level]->protobox_[iter->first][childiter] = iter->second[childiter];

                        }

                    }

                }

            }

            std::unordered_map<long long, std::array<std::unordered_map<long long, std::vector<int>>, 8>>().swap(protobox);

            }

            #pragma omp barrier

            old_cocentered_morton_box = -1;

            #pragma omp for
            for (long long iter = 0; iter < levels_[level-1]->split_points_relconesegments_[comm_rank_+1]-levels_[level-1]->split_points_relconesegments_[comm_rank_]; iter++) {

                const long long cocentered_morton_box = levels_[level-1]->relconesegment2cocenteredmortonboxid_[iter];
                const long long nonrel_conesegment = levels_[level-1]->relconesegment2nonrelconesegment_[iter];

                if (cocentered_morton_box != old_cocentered_morton_box) {

                    morton_children = GetMortonRelChildren(level-1, cocentered_morton_box);

                    old_cocentered_morton_box = cocentered_morton_box;

                }

                for (int childiter = 0; childiter < morton_children.size(); childiter++) {

                    const long long morton_child = morton_children[childiter];
                    const long long child_pos = morton_child % 8;

                    const std::array<long long, 2> data_child = levels_[level]->mortonbox2discretizationpoints_2_[morton_child];

                    if (USE_ADAPTIVITY && (data_child[1] <= MAX_ELEMS_LEAF)) {

                        continue;

                    }

                    const std::unordered_map<long long, std::vector<int>> & allnonrelconesegments = levels_[level]->protobox_.at(nonrel_conesegment)[child_pos];

                    for (auto i = allnonrelconesegments.begin(); i != allnonrelconesegments.end(); i++) {

                        mortonboxnonrelconesegment_thread[morton_child].insert(i->first);

                    }

                }

            }

            #pragma omp critical
            {

            for (auto iter = mortonboxnonrelconesegment_thread.begin(); iter != mortonboxnonrelconesegment_thread.end(); iter++) {

                mortonboxnonrelconesegment_all[iter->first].insert(iter->second.begin(), iter->second.end());

            }

            std::unordered_map<long long, std::unordered_set<long long>>().swap(mortonboxnonrelconesegment_thread);

            }

            }

            for (auto iter = mortonboxnonrelconesegment_all.begin(); iter != mortonboxnonrelconesegment_all.end(); iter++) {

                std::vector<long long> tmp(iter->second.begin(), iter->second.end());
                std::sort(tmp.begin(), tmp.end());
                mortonboxnonrelconesegment[iter->first] = std::move(tmp);

            }

        }

        void GetRelevantConeSegmentsDueToParents(const int level, 
            std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment,
            std::unordered_map<long long, std::unordered_map<long long, double>> & min_radius_conesegment,
            std::unordered_map<long long, std::unordered_map<long long, double>> & max_radius_conesegment) 
        {

            std::unordered_map<long long, std::unordered_set<long long>> mortonboxnonrelconesegment_all;

            #pragma omp parallel
            {

            long long old_cocentered_morton_box = -1;
            std::vector<long long> morton_children;
            std::array<double, 3> center;

            std::unordered_map<long long, std::unordered_set<long long>> mortonboxnonrelconesegment_thread;
            std::unordered_map<long long, std::unordered_map<long long, double>> tmp_min_radius_conesegment;
            std::unordered_map<long long, std::unordered_map<long long, double>> tmp_max_radius_conesegment;

            #pragma omp for
            for (long long iter = 0; iter < levels_[level-1]->split_points_relconesegments_[comm_rank_+1]-levels_[level-1]->split_points_relconesegments_[comm_rank_]; iter++) {

                const long long cocentered_morton_box = levels_[level-1]->relconesegment2cocenteredmortonboxid_[iter];
                const long long nonrel_conesegment = levels_[level-1]->relconesegment2nonrelconesegment_[iter];
                const double radius1 = levels_[level-1]->relconesegment2minradius_[iter];
                const double radius2 = levels_[level-1]->relconesegment2maxradius_[iter];

                Eigen::MatrixXd interp_points = levels_[level-1]->interpolation_points_[nonrel_conesegment];

                if (cocentered_morton_box != old_cocentered_morton_box) {

                    morton_children = GetMortonRelChildren(level-1, cocentered_morton_box);

                    if (morton_children.size() == 0) {
                        throw std::logic_error("Every cocentered morton box must have children in get relevant cone segments function.");
                    }

                    old_cocentered_morton_box = cocentered_morton_box;

                    levels_[level-1]->GetBoxCenter(cocentered_morton_box, center[0], center[1], center[2]);

                }

                for (int childiter = 0; childiter < morton_children.size(); childiter++) {

                    const long long morton_child = morton_children[childiter];

                    std::array<double, 3> center_child;
                    levels_[level]->GetBoxCenter(morton_child, center_child[0], center_child[1], center_child[2]);

                    const std::array<long long, 2> data_child = levels_[level]->mortonbox2discretizationpoints_2_[morton_child];

                    if (USE_ADAPTIVITY && (data_child[1] <= MAX_ELEMS_LEAF)) {

                        continue;

                    }

                    for (int interpiter = 0; interpiter < PT_II*PS_II; interpiter++) {

                        const double x = center[0] + radius1 * interp_points(0, interpiter);
                        const double y = center[1] + radius1 * interp_points(1, interpiter);
                        const double z = center[2] + radius1 * interp_points(2, interpiter);

                        const long long locnonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_child, x, y, z);
                        
                        mortonboxnonrelconesegment_thread[morton_child].insert(locnonrel_conesegment);

                        const double diff = std::sqrt((x-center_child[0])*(x-center_child[0]) + (y-center_child[1])*(y-center_child[1]) + (z-center_child[2])*(z-center_child[2]));

                        if ((tmp_min_radius_conesegment.count(morton_child) == 0) || (tmp_min_radius_conesegment.at(morton_child).count(locnonrel_conesegment) == 0)) {

                            tmp_min_radius_conesegment[morton_child][locnonrel_conesegment] = diff - 1E-10;

                        } else {

                            tmp_min_radius_conesegment[morton_child][locnonrel_conesegment] = std::min({tmp_min_radius_conesegment[morton_child][locnonrel_conesegment], diff - 1E-10});

                        }

                        if ((tmp_max_radius_conesegment.count(morton_child) == 0) || (tmp_max_radius_conesegment.at(morton_child).count(locnonrel_conesegment) == 0)) {

                            tmp_max_radius_conesegment[morton_child][locnonrel_conesegment] = diff + 1E-10;

                        } else {

                            tmp_max_radius_conesegment[morton_child][locnonrel_conesegment] = std::max({tmp_max_radius_conesegment[morton_child][locnonrel_conesegment], diff + 1E-10});

                        }

                    }

                    for (int interpiter = 0; interpiter < PT_II*PS_II; interpiter++) {

                        const double x = center[0] + radius2 * interp_points(0, interpiter);
                        const double y = center[1] + radius2 * interp_points(1, interpiter);
                        const double z = center[2] + radius2 * interp_points(2, interpiter);

                        const long long locnonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_child, x, y, z);
                        
                        mortonboxnonrelconesegment_thread[morton_child].insert(locnonrel_conesegment);

                        const double diff = std::sqrt((x-center_child[0])*(x-center_child[0]) + (y-center_child[1])*(y-center_child[1]) + (z-center_child[2])*(z-center_child[2]));

                        if ((tmp_min_radius_conesegment.count(morton_child) == 0) || (tmp_min_radius_conesegment.at(morton_child).count(locnonrel_conesegment) == 0)) {

                            tmp_min_radius_conesegment[morton_child][locnonrel_conesegment] = diff - 1E-10;

                        } else {

                            tmp_min_radius_conesegment[morton_child][locnonrel_conesegment] = std::min({tmp_min_radius_conesegment[morton_child][locnonrel_conesegment], diff - 1E-10});

                        }

                        if ((tmp_max_radius_conesegment.count(morton_child) == 0) || (tmp_max_radius_conesegment.at(morton_child).count(locnonrel_conesegment) == 0)) {

                            tmp_max_radius_conesegment[morton_child][locnonrel_conesegment] = diff + 1E-10;

                        } else {

                            tmp_max_radius_conesegment[morton_child][locnonrel_conesegment] = std::max({tmp_max_radius_conesegment[morton_child][locnonrel_conesegment], diff + 1E-10});

                        }

                    }

                }

            }

            #pragma omp critical
            {

            for (auto iter = mortonboxnonrelconesegment_thread.begin(); iter != mortonboxnonrelconesegment_thread.end(); iter++) {

                mortonboxnonrelconesegment_all[iter->first].insert(iter->second.begin(), iter->second.end());

            }

            std::unordered_map<long long, std::unordered_set<long long>>().swap(mortonboxnonrelconesegment_thread);

            for (const auto & [key_1, map_1] : tmp_min_radius_conesegment) {
                for (const auto & [key_2, value] : map_1) {

                    if ((min_radius_conesegment.count(key_1) == 0) || (min_radius_conesegment[key_1].count(key_2) == 0)) {

                        min_radius_conesegment[key_1][key_2] = value;

                    } else {

                        min_radius_conesegment[key_1][key_2] = std::min({min_radius_conesegment[key_1][key_2], value});

                    }

                }
            }

            for (const auto & [key_1, map_1] : tmp_max_radius_conesegment) {
                for (const auto & [key_2, value] : map_1) {

                    if ((max_radius_conesegment.count(key_1) == 0) || (max_radius_conesegment[key_1].count(key_2) == 0)) {

                        max_radius_conesegment[key_1][key_2] = value;

                    } else {

                        max_radius_conesegment[key_1][key_2] = std::max({max_radius_conesegment[key_1][key_2], value});

                    }

                }
            }

            std::unordered_map<long long, std::unordered_map<long long, double>>().swap(tmp_min_radius_conesegment);
            std::unordered_map<long long, std::unordered_map<long long, double>>().swap(tmp_max_radius_conesegment);

            }

            }

            for (auto iter = mortonboxnonrelconesegment_all.begin(); iter != mortonboxnonrelconesegment_all.end(); iter++) {

                std::vector<long long> tmp(iter->second.begin(), iter->second.end());
                std::sort(tmp.begin(), tmp.end());
                mortonboxnonrelconesegment[iter->first] = std::move(tmp);

            }          

        }
        
        std::unordered_map<long long, std::vector<long long>> MergeRelevantConeSegments(
            const std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment1,
            const std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment2) const 
        {

            std::unordered_map<long long, std::vector<long long>> merged_map = mortonboxnonrelconesegment1;

            for (auto iter = mortonboxnonrelconesegment2.begin(); iter != mortonboxnonrelconesegment2.end(); iter++) {

                if (merged_map.count(iter->first) == 0) 
                    merged_map.insert({iter->first, std::vector<long long>{}});

            }

            for (auto iter = mortonboxnonrelconesegment2.begin(); iter != mortonboxnonrelconesegment2.end(); iter++) {

                std::vector<long long> & vec1 = merged_map.at(iter->first);

                std::vector<long long> vec;
                vec.reserve(iter->second.size() + vec1.size());
                std::merge(vec1.begin(), vec1.end(), iter->second.begin(), iter->second.end(), std::back_inserter(vec));

                auto last = std::unique(std::begin(vec), std::end(vec));
                vec.erase(last, std::end(vec));
                vec1 = std::move(vec);

            }

            return merged_map;

        }

        void LoadBalanceRelevantConeSegments(const int level,
            const std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment) 
        {

            std::map<long long, std::set<long long>> allmortonboxnonrelconesegment;

            const int nboxes = mortonboxnonrelconesegment.size();

            std::vector<long long> nconesperbox(nboxes);
            std::vector<long long> mortonboxids(nboxes);
            std::vector<long long> nonrelconesegments;

            int boxcounter = 0;
            int size_nonrelconesegments = 0;

            for (const auto & [key, value] : mortonboxnonrelconesegment) {

                nconesperbox[boxcounter] = value.size();
                mortonboxids[boxcounter] = key;
                nonrelconesegments.insert(nonrelconesegments.end(), value.begin(), value.end());

                boxcounter++;
                size_nonrelconesegments += value.size();

            }

            std::vector<int> nboxes_all;
            std::vector<int> size_nonrelconesegments_all;

            if (comm_rank_ == 0) {

                nboxes_all.resize(comm_size_);
                size_nonrelconesegments_all.resize(comm_size_);

            }

            MPI_Gather(&nboxes, 1, MPI_INT, &nboxes_all[0], 1, MPI_INT, 0, mpi_comm_);
            MPI_Gather(&size_nonrelconesegments, 1, MPI_INT, &size_nonrelconesegments_all[0], 1, MPI_INT, 0, mpi_comm_);

            std::vector<long long> nconesperbox_all;
            std::vector<long long> mortonboxids_all;
            std::vector<long long> nonrelconesegments_all;
            std::vector<int> displs;
            std::vector<int> displs_2;

            if (comm_rank_ == 0) {     

                long long nboxes_total = 0;
                long long nonrelconesegments_total = 0;

                for (int i = 0; i < comm_size_; i++) {

                    nboxes_total += nboxes_all[i];
                    nonrelconesegments_total += size_nonrelconesegments_all[i];

                }

                nconesperbox_all.resize(nboxes_total);
                mortonboxids_all.resize(nboxes_total);
                nonrelconesegments_all.resize(nonrelconesegments_total);
                displs.resize(comm_size_);
                displs_2.resize(comm_size_);

                displs[0] = 0;
                displs_2[0] = 0;

                for (int i = 1; i < comm_size_; i++) {

                    displs[i] = nboxes_all[i-1] + displs[i-1];
                    displs_2[i] = size_nonrelconesegments_all[i-1] + displs_2[i-1];

                }

            }

            MPI_Barrier(mpi_comm_);
            
            MPI_Gatherv(&nconesperbox[0], nboxes, MPI_LONG_LONG, &nconesperbox_all[0], &nboxes_all[0], &displs[0], MPI_LONG_LONG, 0, mpi_comm_);
            MPI_Gatherv(&mortonboxids[0], nboxes, MPI_LONG_LONG, &mortonboxids_all[0], &nboxes_all[0], &displs[0], MPI_LONG_LONG, 0, mpi_comm_);
            MPI_Gatherv(&nonrelconesegments[0], size_nonrelconesegments, MPI_LONG_LONG, &nonrelconesegments_all[0], &size_nonrelconesegments_all[0], &displs_2[0], MPI_LONG_LONG, 0, mpi_comm_);

            std::vector<long long>().swap(nconesperbox);
            std::vector<long long>().swap(mortonboxids);
            std::vector<long long>().swap(nonrelconesegments);

            if (comm_rank_ == 0) {

                long long idx = 0;

                for (long long i = 0; i < mortonboxids_all.size(); i++) {

                    long long morton_box = mortonboxids_all[i];

                    long long size = nconesperbox_all[i];

                    for (long long j = 0; j < size; j++) {

                        allmortonboxnonrelconesegment[morton_box].insert(nonrelconesegments_all[idx]);
                        idx++;

                    }                    

                }

            }

            std::vector<long long>().swap(nconesperbox_all);
            std::vector<long long>().swap(mortonboxids_all);
            std::vector<long long>().swap(nonrelconesegments_all);

            long long nallconesegments = 0;

            std::vector<long long> relcone2nonrelcone;
            std::vector<long long> relcone2cocenteredmorton;

            if (comm_rank_ == 0) {

                for (const auto & [key, value] : allmortonboxnonrelconesegment) {

                    nallconesegments += value.size();

                }

                relcone2nonrelcone.resize(nallconesegments);
                relcone2cocenteredmorton.resize(nallconesegments);

                long long counter = 0;

                for (const auto & [morton_box, cone_set] : allmortonboxnonrelconesegment) {
                    for (const auto & cone : cone_set) {

                        relcone2nonrelcone[counter] = cone;
                        relcone2cocenteredmorton[counter] = morton_box;
                        counter++;

                    }
                }

            }

            allmortonboxnonrelconesegment.clear();

            MPI_Bcast(&nallconesegments, 1, MPI_LONG_LONG, 0, mpi_comm_);

            if (comm_rank_ != 0) {

                relcone2nonrelcone.resize(nallconesegments);
                relcone2cocenteredmorton.resize(nallconesegments);

            }

            MPI_Bcast(&relcone2nonrelcone[0], nallconesegments, MPI_LONG_LONG, 0, mpi_comm_);
            MPI_Bcast(&relcone2cocenteredmorton[0], nallconesegments, MPI_LONG_LONG, 0, mpi_comm_);

            levels_[level]->split_points_relconesegments_.resize(comm_size_+1);

            const long long conesegments_per_rank = nallconesegments / comm_size_;
            long long remaining_conesegments = nallconesegments % comm_size_;

            levels_[level]->split_points_relconesegments_[0] = 0;
            levels_[level]->split_points_relconesegments_[comm_size_] = nallconesegments;

            for (int i = 1; i < comm_size_; i++) {

                levels_[level]->split_points_relconesegments_[i] = levels_[level]->split_points_relconesegments_[i-1] + conesegments_per_rank;

                if (remaining_conesegments > 0) {
                    levels_[level]->split_points_relconesegments_[i]++;
                    remaining_conesegments--;
                }

            }

            const long long conesegments_low = levels_[level]->split_points_relconesegments_[comm_rank_];
            const long long conesegments_up = levels_[level]->split_points_relconesegments_[comm_rank_+1];

            levels_[level]->relconesegment2nonrelconesegment_ = std::vector<long long>(conesegments_up-conesegments_low);
            levels_[level]->relconesegment2cocenteredmortonboxid_ = std::vector<long long>(conesegments_up-conesegments_low);

            for (long long i = 0; i < conesegments_up-conesegments_low; i++) {

                levels_[level]->relconesegment2nonrelconesegment_[i] = relcone2nonrelcone[i + conesegments_low];
                levels_[level]->relconesegment2cocenteredmortonboxid_[i] = relcone2cocenteredmorton[i + conesegments_low];
                levels_[level]->mortonboxnonrelcone2relcone_[relcone2cocenteredmorton[i + conesegments_low]][relcone2nonrelcone[i + conesegments_low]] = i + conesegments_low;

            }
            
            for (const auto & [morton_box, cones_vec] : mortonboxnonrelconesegment) {

                for (const auto & cone : cones_vec) {

                    auto box_begin = std::lower_bound(relcone2cocenteredmorton.begin(), relcone2cocenteredmorton.end(), morton_box);
                    auto box_end = std::upper_bound(relcone2cocenteredmorton.begin(), relcone2cocenteredmorton.end(), morton_box);
                    auto cones_begin = relcone2nonrelcone.begin() + std::distance(relcone2cocenteredmorton.begin(), box_begin);
                    auto cones_end = relcone2nonrelcone.begin() + std::distance(relcone2cocenteredmorton.begin(), box_end);
                    auto pos = std::lower_bound(cones_begin, cones_end, cone);

                    if (pos == cones_end) {

                        std::cout << std::endl << "RANK = " << comm_rank_ << ", mortonbox, nonrelcone = " << morton_box << ", " << cone << std::endl;
                        throw std::logic_error("The cone has to be in this box");
                    
                    }

                    levels_[level]->mortonboxnonrelcone2relcone_2_[morton_box][cone] = std::distance(relcone2nonrelcone.begin(), pos);

                }

            }            

            std::vector<long long>().swap(relcone2nonrelcone);
            std::vector<long long>().swap(relcone2cocenteredmorton);

        }

        std::unordered_map<long long, std::unordered_map<long long, double>> MergeRelevantConeSegmentsMinRadius(
            const std::unordered_map<long long, std::unordered_map<long long, double>> & min_radius1,
            const std::unordered_map<long long, std::unordered_map<long long, double>> & min_radius2) const 
        {

            std::unordered_map<long long, std::unordered_map<long long, double>> merged_map = min_radius1;

            for (const auto & [morton_box, map_1] : min_radius2) {

                for (const auto & [cone_segment, radius] : map_1) {

                    if (merged_map.count(morton_box) == 0 || merged_map.at(morton_box).count(cone_segment) == 0) {

                        merged_map[morton_box][cone_segment] = radius;

                    } else {

                        const double original_value = merged_map.at(morton_box).at(cone_segment);

                        merged_map[morton_box][cone_segment] = std::min({original_value, radius});

                    }

                }

            }

            return merged_map;

        }

        std::unordered_map<long long, std::unordered_map<long long, double>> MergeRelevantConeSegmentsMaxRadius(
            const std::unordered_map<long long, std::unordered_map<long long, double>> & max_radius1,
            const std::unordered_map<long long, std::unordered_map<long long, double>> & max_radius2) const 
        {

            std::unordered_map<long long, std::unordered_map<long long, double>> merged_map = max_radius1;

            for (const auto & [morton_box, map_1] : max_radius2) {

                for (const auto & [cone_segment, radius] : map_1) {

                    if (merged_map.count(morton_box) == 0 || merged_map.at(morton_box).count(cone_segment) == 0) {

                        merged_map[morton_box][cone_segment] = radius;

                    } else {

                        const double original_value = merged_map.at(morton_box).at(cone_segment);

                        merged_map[morton_box][cone_segment] = std::max({original_value, radius});

                    }

                }

            }

            return merged_map;

        }

        void LoadBalanceRelevantConeSegmentsMinRadius(const int level,
            const std::unordered_map<long long, std::vector<long long>> & mortonboxnonrelconesegment,
            const std::unordered_map<long long, std::unordered_map<long long, double>> & min_radius_conesegments,
            const std::unordered_map<long long, std::unordered_map<long long, double>> & max_radius_conesegments) 
        {

            std::map<long long, std::set<long long>> allmortonboxnonrelconesegment;
            std::unordered_map<long long, std::unordered_map<long long, double>> all_min_radius_conesegments;
            std::unordered_map<long long, std::unordered_map<long long, double>> all_max_radius_conesegments;
            
            const int nboxes = mortonboxnonrelconesegment.size();

            std::vector<long long> nconesperbox(nboxes);
            std::vector<long long> mortonboxids(nboxes);
            std::vector<long long> nonrelconesegments;
            std::vector<double> min_radius;
            std::vector<double> max_radius;

            int boxcounter = 0;
            int size_nonrelconesegments = 0;

            for (const auto & [key, value] : mortonboxnonrelconesegment) {

                nconesperbox[boxcounter] = value.size();
                mortonboxids[boxcounter] = key;
                nonrelconesegments.insert(nonrelconesegments.end(), value.begin(), value.end());

                for (long long i = 0; i < value.size(); i++) {

                    min_radius.push_back(min_radius_conesegments.at(key).at(value[i]));
                    max_radius.push_back(max_radius_conesegments.at(key).at(value[i]));

                }

                boxcounter++;
                size_nonrelconesegments += value.size();

            }

            std::vector<int> nboxes_all;
            std::vector<int> size_nonrelconesegments_all;

            if (comm_rank_ == 0) {

                nboxes_all.resize(comm_size_);
                size_nonrelconesegments_all.resize(comm_size_);

            }

            MPI_Gather(&nboxes, 1, MPI_INT, &nboxes_all[0], 1, MPI_INT, 0, mpi_comm_);
            MPI_Gather(&size_nonrelconesegments, 1, MPI_INT, &size_nonrelconesegments_all[0], 1, MPI_INT, 0, mpi_comm_);

            std::vector<long long> nconesperbox_all;
            std::vector<long long> mortonboxids_all;
            std::vector<long long> nonrelconesegments_all;
            std::vector<double> min_radius_all;
            std::vector<double> max_radius_all;
            std::vector<int> displs;
            std::vector<int> displs_2;

            if (comm_rank_ == 0) {     

                long long nboxes_total = 0;
                long long nonrelconesegments_total = 0;

                for (int i = 0; i < comm_size_; i++) {

                    nboxes_total += nboxes_all[i];
                    nonrelconesegments_total += size_nonrelconesegments_all[i];

                }

                nconesperbox_all.resize(nboxes_total);
                mortonboxids_all.resize(nboxes_total);
                nonrelconesegments_all.resize(nonrelconesegments_total);
                min_radius_all.resize(nonrelconesegments_total);
                max_radius_all.resize(nonrelconesegments_total);
                displs.resize(comm_size_);
                displs_2.resize(comm_size_);
                
                displs[0] = 0;
                displs_2[0] = 0;

                for (int i = 1; i < comm_size_; i++) {

                    displs[i] = nboxes_all[i-1] + displs[i-1];
                    displs_2[i] = size_nonrelconesegments_all[i-1] + displs_2[i-1];
                    
                }

            }

            MPI_Barrier(mpi_comm_);
            
            MPI_Gatherv(&nconesperbox[0], nboxes, MPI_LONG_LONG, &nconesperbox_all[0], &nboxes_all[0], &displs[0], MPI_LONG_LONG, 0, mpi_comm_);
            MPI_Gatherv(&mortonboxids[0], nboxes, MPI_LONG_LONG, &mortonboxids_all[0], &nboxes_all[0], &displs[0], MPI_LONG_LONG, 0, mpi_comm_);
            MPI_Gatherv(&nonrelconesegments[0], size_nonrelconesegments, MPI_LONG_LONG, &nonrelconesegments_all[0], &size_nonrelconesegments_all[0], &displs_2[0], MPI_LONG_LONG, 0, mpi_comm_);
            MPI_Gatherv(&min_radius[0], size_nonrelconesegments, MPI_DOUBLE, &min_radius_all[0], &size_nonrelconesegments_all[0], &displs_2[0], MPI_DOUBLE, 0, mpi_comm_);
            MPI_Gatherv(&max_radius[0], size_nonrelconesegments, MPI_DOUBLE, &max_radius_all[0], &size_nonrelconesegments_all[0], &displs_2[0], MPI_DOUBLE, 0, mpi_comm_);

            std::vector<long long>().swap(nconesperbox);
            std::vector<long long>().swap(mortonboxids);
            std::vector<long long>().swap(nonrelconesegments);
            std::vector<double>().swap(min_radius);
            std::vector<double>().swap(max_radius);

            if (comm_rank_ == 0) {

                long long idx = 0;

                for (long long i = 0; i < mortonboxids_all.size(); i++) {

                    const long long morton_box = mortonboxids_all[i];

                    const long long size = nconesperbox_all[i];

                    for (long long j = 0; j < size; j++) {

                        const long long cone_segment = nonrelconesegments_all[idx];

                        allmortonboxnonrelconesegment[morton_box].insert(cone_segment);

                        const double radius = min_radius_all[idx];

                        if (all_min_radius_conesegments.count(morton_box) == 0 || all_min_radius_conesegments[morton_box].count(cone_segment) == 0) {

                            all_min_radius_conesegments[morton_box][cone_segment] = radius;

                        } else {

                            const double original_value = all_min_radius_conesegments[morton_box][cone_segment];

                            all_min_radius_conesegments[morton_box][cone_segment] = std::min({original_value, radius});

                        }    

                        const double radius2 = max_radius_all[idx];

                        if (all_max_radius_conesegments.count(morton_box) == 0 || all_max_radius_conesegments[morton_box].count(cone_segment) == 0) {

                            all_max_radius_conesegments[morton_box][cone_segment] = radius2;

                        } else {

                            const double original_value = all_max_radius_conesegments[morton_box][cone_segment];

                            all_max_radius_conesegments[morton_box][cone_segment] = std::max({original_value, radius2});

                        }                    
                        
                        idx++;

                    }                       

                }

            }
        
            std::vector<long long>().swap(nconesperbox_all);
            std::vector<long long>().swap(mortonboxids_all);
            std::vector<long long>().swap(nonrelconesegments_all);
            std::vector<double>().swap(min_radius_all);
            std::vector<double>().swap(max_radius_all);

            long long nallconesegments = 0;

            std::vector<long long> relcone2nonrelcone;
            std::vector<long long> relcone2cocenteredmorton;
            std::vector<double> relcone2minradius;
            std::vector<double> relcone2maxradius;

            if (comm_rank_ == 0) {

                for (const auto & [key, value] : allmortonboxnonrelconesegment) {

                    nallconesegments += value.size();

                }

                relcone2nonrelcone.resize(nallconesegments);
                relcone2cocenteredmorton.resize(nallconesegments);
                relcone2minradius.resize(nallconesegments);
                relcone2maxradius.resize(nallconesegments);

                long long counter = 0;

                for (const auto & [morton_box, cone_set] : allmortonboxnonrelconesegment) {
                    for (const auto & cone : cone_set) {

                        relcone2nonrelcone[counter] = cone;
                        relcone2cocenteredmorton[counter] = morton_box;
                        relcone2minradius[counter] = all_min_radius_conesegments.at(morton_box).at(cone);
                        relcone2maxradius[counter] = all_max_radius_conesegments.at(morton_box).at(cone);
                        
                        counter++;

                    }
                }

            }

            allmortonboxnonrelconesegment.clear();
            all_min_radius_conesegments.clear();
            all_max_radius_conesegments.clear();

            MPI_Bcast(&nallconesegments, 1, MPI_LONG_LONG, 0, mpi_comm_);

            if (comm_rank_ != 0) {

                relcone2nonrelcone.resize(nallconesegments);
                relcone2cocenteredmorton.resize(nallconesegments);
                relcone2minradius.resize(nallconesegments);
                relcone2maxradius.resize(nallconesegments);

            }

            MPI_Bcast(&relcone2nonrelcone[0], nallconesegments, MPI_LONG_LONG, 0, mpi_comm_);
            MPI_Bcast(&relcone2cocenteredmorton[0], nallconesegments, MPI_LONG_LONG, 0, mpi_comm_);
            MPI_Bcast(&relcone2minradius[0], nallconesegments, MPI_DOUBLE, 0, mpi_comm_);
            MPI_Bcast(&relcone2maxradius[0], nallconesegments, MPI_DOUBLE, 0, mpi_comm_);

            levels_[level]->split_points_relconesegments_.resize(comm_size_+1);

            const long long conesegments_per_rank = nallconesegments / comm_size_;
            long long remaining_conesegments = nallconesegments % comm_size_;

            levels_[level]->split_points_relconesegments_[0] = 0;
            levels_[level]->split_points_relconesegments_[comm_size_] = nallconesegments;

            for (int i = 1; i < comm_size_; i++) {

                levels_[level]->split_points_relconesegments_[i] = levels_[level]->split_points_relconesegments_[i-1] + conesegments_per_rank;

                if (remaining_conesegments > 0) {
                    levels_[level]->split_points_relconesegments_[i]++;
                    remaining_conesegments--;
                }

            }

            const long long conesegments_low = levels_[level]->split_points_relconesegments_[comm_rank_];
            const long long conesegments_up = levels_[level]->split_points_relconesegments_[comm_rank_+1];

            levels_[level]->relconesegment2nonrelconesegment_ = std::vector<long long>(conesegments_up-conesegments_low);
            levels_[level]->relconesegment2cocenteredmortonboxid_ = std::vector<long long>(conesegments_up-conesegments_low);
            levels_[level]->relconesegment2minradius_ = std::vector<double>(conesegments_up-conesegments_low);
            levels_[level]->relconesegment2maxradius_ = std::vector<double>(conesegments_up-conesegments_low);
            
            for (long long i = 0; i < conesegments_up-conesegments_low; i++) {

                levels_[level]->relconesegment2nonrelconesegment_[i] = relcone2nonrelcone[i + conesegments_low];
                levels_[level]->relconesegment2cocenteredmortonboxid_[i] = relcone2cocenteredmorton[i + conesegments_low];
                levels_[level]->mortonboxnonrelcone2relcone_[relcone2cocenteredmorton[i + conesegments_low]][relcone2nonrelcone[i + conesegments_low]] = i + conesegments_low;
                levels_[level]->relconesegment2minradius_[i] = relcone2minradius[i + conesegments_low];
                levels_[level]->relconesegment2maxradius_[i] = relcone2maxradius[i + conesegments_low];

            }

            for (const auto & [morton_box, cones_vec] : mortonboxnonrelconesegment) {

                for (const auto & cone : cones_vec) {

                    auto box_begin = std::lower_bound(relcone2cocenteredmorton.begin(), relcone2cocenteredmorton.end(), morton_box);
                    auto box_end = std::upper_bound(relcone2cocenteredmorton.begin(), relcone2cocenteredmorton.end(), morton_box);
                    auto cones_begin = relcone2nonrelcone.begin() + std::distance(relcone2cocenteredmorton.begin(), box_begin);
                    auto cones_end = relcone2nonrelcone.begin() + std::distance(relcone2cocenteredmorton.begin(), box_end);
                    auto pos = std::lower_bound(cones_begin, cones_end, cone);

                    if (pos == cones_end) {

                        std::cout << std::endl << "RANK = " << comm_rank_ << ", mortonbox, nonrelcone = " << morton_box << ", " << cone << std::endl;
                        throw std::logic_error("The cone has to be in this box");
                    
                    }

                    long long relconesegment = std::distance(relcone2nonrelcone.begin(), pos);

                    levels_[level]->mortonboxnonrelcone2relcone_2_[morton_box][cone] = relconesegment;

                    levels_[level]->mortonboxnonrelcone2minradius_[relconesegment] = relcone2minradius[relconesegment];
                    levels_[level]->mortonboxnonrelcone2maxradius_[relconesegment] = relcone2maxradius[relconesegment];

                }

            }              

            std::vector<long long>().swap(relcone2nonrelcone);
            std::vector<long long>().swap(relcone2cocenteredmorton);
            std::vector<double>().swap(relcone2minradius);
            std::vector<double>().swap(relcone2maxradius);

        }

        void SetRelevantChildrenOfCocenteredBoxes(int level)
        {

            if (level < nlevels_ - 1) {

                std::vector<long long> morton_box_loc;
                std::vector<long long> position_loc;
                std::vector<long long> size_loc;

                for (const auto & [morton_box, values] : levels_[level+1]->mortonbox2discretizationpoints_) {

                    morton_box_loc.push_back(morton_box);
                    position_loc.push_back(values[0]);
                    size_loc.push_back(values[1]);

                }

                int total_size_loc = morton_box_loc.size();

                std::vector<int> recv_counts(comm_size_);
                MPI_Allgather(&total_size_loc, 1, MPI_INT, &recv_counts[0], 1, MPI_INT, mpi_comm_);

                std::vector<int> displs(comm_size_, 0);
                long long total_size = 0;

                for (int i = 0; i < comm_size_; i++) {

                    total_size += recv_counts[i];

                    if (i != 0) {
                        displs[i] = displs[i-1] + recv_counts[i-1];
                    }

                }

                std::vector<long long> morton_box_all(total_size);
                std::vector<long long> position_all(total_size);
                std::vector<long long> size_all(total_size);

                MPI_Allgatherv(&morton_box_loc[0], total_size_loc, MPI_LONG_LONG, &morton_box_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
                MPI_Allgatherv(&position_loc[0], total_size_loc, MPI_LONG_LONG, &position_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
                MPI_Allgatherv(&size_loc[0], total_size_loc, MPI_LONG_LONG, &size_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);

                std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_all;

                for (long long i = 0; i < total_size; i++) {

                    const long long morton_box = morton_box_all[i];
                    const long long position = position_all[i];
                    const long long size = size_all[i];

                    if (mortonbox2discretizationpoints_all.count(morton_box) == 0) {

                        mortonbox2discretizationpoints_all[morton_box] = {position, size};

                    } else {

                        const std::array<long long, 2> actual_value = mortonbox2discretizationpoints_all[morton_box];

                        mortonbox2discretizationpoints_all[morton_box] = {std::min({actual_value[0], position}), actual_value[1] + size};

                    }

                }

                std::vector<long long>().swap(morton_box_loc);
                std::vector<long long>().swap(position_loc);
                std::vector<long long>().swap(size_loc);

                std::vector<long long>().swap(morton_box_all);
                std::vector<long long>().swap(position_all);
                std::vector<long long>().swap(size_all);

                for (const auto & [morton_box, cones] : levels_[level]->mortonboxnonrelcone2relcone_) {

                    for (long long childiter = morton_box*8; childiter < (morton_box+1)*8; childiter++) {

                        if (mortonbox2discretizationpoints_all.count(childiter) != 0 && levels_[level+1]->mortonidofrelboxes_2_.count(childiter) == 0) {
                           
                            levels_[level+1]->mortonidofrelboxes_2_.insert(childiter);
                            levels_[level+1]->mortonbox2discretizationpoints_2_[childiter] = mortonbox2discretizationpoints_all[childiter];
                        
                        }

                    }

                }

            }

        } 

        void SetPropagationAndInterpolationRequiredConeSegments(const int level, 
            const std::unordered_map<long long, std::vector<long long>> & interpolation_conesegments,
            const std::unordered_map<long long, std::vector<long long>> & propagation_conesegments) 
        {

            levels_[level]->interpolationrequiredrelconesegments_.clear();
            levels_[level]->propagationrequiredrelconesegments_.clear();

            for (const auto & [morton_box, cone_segments] : interpolation_conesegments) {

                for (long long csiter = 0; csiter < cone_segments.size(); csiter++) {

                    const long long relconesegment = levels_[level]->MortonboxNonrelconesegment2Relconesegment(morton_box, cone_segments[csiter]);

                    levels_[level]->interpolationrequiredrelconesegments_.insert(relconesegment);

                }

            }

            for (const auto & [morton_box, cone_segments] : propagation_conesegments) {

                for (long long csiter = 0; csiter < cone_segments.size(); csiter++) {

                    const long long relconesegment = levels_[level]->MortonboxNonrelconesegment2Relconesegment(morton_box, cone_segments[csiter]);

                    levels_[level]->propagationrequiredrelconesegments_.insert(relconesegment);

                }

            }

        }

        void SetRelevantConeSegmentsNotInRank(int level) {

            levels_[level]->propagation_not_in_rank_ = std::vector<std::vector<long long>>(comm_size_);
            levels_[level]->interpolation_not_in_rank_ = std::vector<std::vector<long long>>(comm_size_);

            for (const auto & cone_segment : levels_[level]->propagationrequiredrelconesegments_) {

                int rank;

                for (int i = 0; i < comm_size_+1; i++) {

                    if ((cone_segment >= levels_[level]->split_points_relconesegments_[i]) && (cone_segment < levels_[level]->split_points_relconesegments_[i+1])) {

                        rank = i;
                        break;

                    }

                }

                if (rank != comm_rank_) {

                    levels_[level]->propagation_not_in_rank_[rank].push_back(cone_segment);

                }

            }

            for (int rank = 0; rank < comm_size_; rank++) {

                std::vector<long long> cone_segments_loc = levels_[level]->propagation_not_in_rank_[rank];
                std::vector<long long> cone_segments_all;

                int size_loc = cone_segments_loc.size();

                std::vector<int> recv_counts(comm_size_);
                MPI_Gather(&size_loc, 1, MPI_INT, &recv_counts[0], 1, MPI_INT, rank, mpi_comm_);

                std::vector<int> displs(comm_size_, 0);
                long long size_all = 0;

                if (comm_rank_ == rank) {

                    size_all = recv_counts[0];

                    for (int i = 1; i < comm_size_; i++) {

                        displs[i] = displs[i-1] + recv_counts[i-1];

                        size_all += recv_counts[i];

                    }
                    
                    cone_segments_all.resize(size_all);

                }

                MPI_Gatherv(&cone_segments_loc[0], size_loc, MPI_LONG_LONG, &cone_segments_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, rank, mpi_comm_);

                if (comm_rank_ == rank) {

                    levels_[level]->recv_counts_propagation_.resize(comm_size_);
                    levels_[level]->displs_propagation_.resize(comm_size_);

                    for (int i = 0; i < comm_size_; i++) {

                        levels_[level]->recv_counts_propagation_[i] = recv_counts[i] * P_;
                        levels_[level]->displs_propagation_[i] = displs[i] * P_;

                    }

                    levels_[level]->cone_segments_propagation_not_in_rank_all_ = cone_segments_all;

                }

            }

            for (const auto & cone_segment : levels_[level]->interpolationrequiredrelconesegments_) {

                int rank;

                for (int i = 0; i < comm_size_+1; i++) {

                    if ((cone_segment >= levels_[level]->split_points_relconesegments_[i]) && (cone_segment < levels_[level]->split_points_relconesegments_[i+1])) {

                        rank = i;
                        break;

                    }

                }

                if (rank != comm_rank_) {

                    levels_[level]->interpolation_not_in_rank_[rank].push_back(cone_segment);

                }

            }

            for (int rank = 0; rank < comm_size_; rank++) {

                std::vector<long long> cone_segments_loc = levels_[level]->interpolation_not_in_rank_[rank];
                std::vector<long long> cone_segments_all;

                int size_loc = cone_segments_loc.size();

                std::vector<int> recv_counts(comm_size_);
                MPI_Gather(&size_loc, 1, MPI_INT, &recv_counts[0], 1, MPI_INT, rank, mpi_comm_);

                std::vector<int> displs(comm_size_, 0);
                long long size_all = 0;

                if (comm_rank_ == rank) {

                    size_all = recv_counts[0];

                    for (int i = 1; i < comm_size_; i++) {

                        displs[i] = displs[i-1] + recv_counts[i-1];

                        size_all += recv_counts[i];

                    }
                    
                    cone_segments_all.resize(size_all);

                }

                MPI_Gatherv(&cone_segments_loc[0], size_loc, MPI_LONG_LONG, &cone_segments_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, rank, mpi_comm_);

                if (comm_rank_ == rank) {

                    levels_[level]->recv_counts_interpolation_.resize(comm_size_);
                    levels_[level]->displs_interpolation_.resize(comm_size_);

                    for (int i = 0; i < comm_size_; i++) {

                        levels_[level]->recv_counts_interpolation_[i] = recv_counts[i] * P_;
                        levels_[level]->displs_interpolation_[i] = displs[i] * P_;

                    }

                    levels_[level]->cone_segments_interpolation_not_in_rank_all_ = cone_segments_all;

                }

            }

        }

        void UpdateBoxesData(int level) 
        {

            std::vector<long long> morton_box_loc;
            std::vector<long long> position_loc;
            std::vector<long long> size_loc;

            for (const auto & [morton_box, values] : levels_[level]->mortonbox2discretizationpoints_) {

                morton_box_loc.push_back(morton_box);
                position_loc.push_back(values[0]);
                size_loc.push_back(values[1]);

            }

            int total_size_loc = morton_box_loc.size();

            std::vector<int> recv_counts(comm_size_);
            MPI_Allgather(&total_size_loc, 1, MPI_INT, &recv_counts[0], 1, MPI_INT, mpi_comm_);

            std::vector<int> displs(comm_size_, 0);
            long long total_size = 0;

            for (int i = 0; i < comm_size_; i++) {

                total_size += recv_counts[i];

                if (i != 0) {
                    displs[i] = displs[i-1] + recv_counts[i-1];
                }

            }

            std::vector<long long> morton_box_all(total_size);
            std::vector<long long> position_all(total_size);
            std::vector<long long> size_all(total_size);

            MPI_Allgatherv(&morton_box_loc[0], total_size_loc, MPI_LONG_LONG, &morton_box_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
            MPI_Allgatherv(&position_loc[0], total_size_loc, MPI_LONG_LONG, &position_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
            MPI_Allgatherv(&size_loc[0], total_size_loc, MPI_LONG_LONG, &size_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);

            std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_all;

            for (long long i = 0; i < total_size; i++) {

                const long long morton_box = morton_box_all[i];
                const long long position = position_all[i];
                const long long size = size_all[i];

                if (mortonbox2discretizationpoints_all.count(morton_box) == 0) {

                    mortonbox2discretizationpoints_all[morton_box] = {position, size};

                } else {

                    const std::array<long long, 2> actual_value = mortonbox2discretizationpoints_all[morton_box];

                    mortonbox2discretizationpoints_all[morton_box] = {std::min({actual_value[0], position}), actual_value[1] + size};

                }

            }

            std::vector<long long>().swap(morton_box_loc);
            std::vector<long long>().swap(position_loc);
            std::vector<long long>().swap(size_loc);

            std::vector<long long>().swap(morton_box_all);
            std::vector<long long>().swap(position_all);
            std::vector<long long>().swap(size_all);
                
            for (const long long & morton_box : levels_[level]->relconesegment2cocenteredmortonboxid_) {
                
                if (levels_[level]->mortonbox2discretizationpoints_2_.count(morton_box) == 0) {

                    levels_[level]->mortonbox2discretizationpoints_2_[morton_box] = mortonbox2discretizationpoints_all[morton_box];
                
                }

            }  

        }

        void InitializeRelevantConeSegments() 
        {

            std::unordered_map<long long, std::vector<long long>> interpolation_conesegments;
            std::unordered_map<long long, std::vector<long long>> propagation_conesegments;

            for (int level = 2; level < nlevels_; level++) {

                GetRelevantConeSegmentsDueToCousinSurface(level, interpolation_conesegments);
                
                if (level > 2) {

                    GetRelevantConeSegmentsDueToParents(level, propagation_conesegments);
                    
                }

                if (interpolation_conesegments.size() + propagation_conesegments.size() == 0) {
                    throw std::logic_error("The total number of cone segments is 0");
                }
                MPI_Barrier(mpi_comm_);                  

                LoadBalanceRelevantConeSegments(level, MergeRelevantConeSegments(interpolation_conesegments, propagation_conesegments));
                MPI_Barrier(mpi_comm_);

                SetRelevantChildrenOfCocenteredBoxes(level);
                MPI_Barrier(mpi_comm_);

                SetPropagationAndInterpolationRequiredConeSegments(level, interpolation_conesegments, propagation_conesegments);
                MPI_Barrier(mpi_comm_);

                SetRelevantConeSegmentsNotInRank(level);
                MPI_Barrier(mpi_comm_);

                UpdateBoxesData(level);
                MPI_Barrier(mpi_comm_);

                interpolation_conesegments.clear();
                propagation_conesegments.clear();
                MPI_Barrier(mpi_comm_);

            }

        }  

        void InitializeRelevantConeSegmentsHO() 
        {

            std::unordered_map<long long, std::vector<long long>> interpolation_conesegments;
            std::unordered_map<long long, std::vector<long long>> propagation_conesegments;

            std::unordered_map<long long, std::unordered_map<long long, double>> min_radius_interpolation_conesegments;
            std::unordered_map<long long, std::unordered_map<long long, double>> min_radius_propagation_conesegments;

            std::unordered_map<long long, std::unordered_map<long long, double>> max_radius_interpolation_conesegments;
            std::unordered_map<long long, std::unordered_map<long long, double>> max_radius_propagation_conesegments;

            for (int level = 2; level < nlevels_; level++) {

                GetRelevantConeSegmentsDueToCousinSurface(level, interpolation_conesegments, min_radius_interpolation_conesegments, max_radius_interpolation_conesegments);

                if (level > 2) {

                    GetRelevantConeSegmentsDueToParents(level, propagation_conesegments, min_radius_propagation_conesegments, max_radius_propagation_conesegments);

                }

                if (interpolation_conesegments.size() + propagation_conesegments.size() == 0) {
                    throw std::logic_error("The total number of cone segments is 0");
                }
                MPI_Barrier(mpi_comm_);
                
                LoadBalanceRelevantConeSegmentsMinRadius(level, MergeRelevantConeSegments(interpolation_conesegments, propagation_conesegments), 
                                                                MergeRelevantConeSegmentsMinRadius(min_radius_interpolation_conesegments, min_radius_propagation_conesegments),
                                                                MergeRelevantConeSegmentsMaxRadius(max_radius_interpolation_conesegments, max_radius_propagation_conesegments));
                MPI_Barrier(mpi_comm_);
                
                LevelPoints(level);
                MPI_Barrier(mpi_comm_);
 
                SetRelevantChildrenOfCocenteredBoxes(level);
                MPI_Barrier(mpi_comm_);

                SetPropagationAndInterpolationRequiredConeSegments(level, interpolation_conesegments, propagation_conesegments);
                MPI_Barrier(mpi_comm_);

                SetRelevantConeSegmentsNotInRank(level);
                MPI_Barrier(mpi_comm_);

                UpdateBoxesData(level);
                MPI_Barrier(mpi_comm_);

                interpolation_conesegments.clear();
                propagation_conesegments.clear();
                min_radius_interpolation_conesegments.clear();
                min_radius_propagation_conesegments.clear();
                MPI_Barrier(mpi_comm_);

            }

        }    

        void InitializeSolutionAndConeSegmentCoefficients() 
        {

            solution_ = std::vector<std::complex<double>>(N_loc_orig_, {0.0, 0.0});
            
            if (USE_HIGH_ORDER) {

                for (int level = 2; level < nlevels_; level++) {

                    long long maxncoeffs = levels_[level]->relconesegment2nonrelconesegment_.size();

                    levels_[level]->indexes_relconesegment_ = std::vector<long long>(maxncoeffs * P_, 0);
                    levels_[level]->num_elems_indexes_relconesegment_ = std::vector<long long>(maxncoeffs);

                }

            }

        }

        void LevelPoints(int level)
        {

            std::unordered_set<long long> unique_nonrel_conesegments(levels_[level]->relconesegment2nonrelconesegment_.begin(), levels_[level]->relconesegment2nonrelconesegment_.end());
            std::vector<long long> unique_nonrel_conesegments_vec(unique_nonrel_conesegments.begin(), unique_nonrel_conesegments.end());

            unique_nonrel_conesegments.clear();

            #pragma omp parallel
            {

            std::unordered_map<long long, Eigen::MatrixXd> interpolation_points_thread;

            #pragma omp for
            for (long long i = 0; i < unique_nonrel_conesegments_vec.size(); i++) {

                const long long nonrel_conesegment = unique_nonrel_conesegments_vec[i];

                long long rr, tt, pp;
                levels_[level]->Nonrelconesegment2threedimensionalconesegment(nonrel_conesegment, rr, tt, pp);

                const double dangle = M_PI / levels_[level]->nconeselevation_;

                const double t_0 = tt * dangle;
                const double t_1 = (tt + 1.0) * dangle;
                const double p_0 = pp * dangle;
                const double p_1 = (pp + 1.0) * dangle;

                const double t_mid = 0.5 * (t_0 + t_1);
                const double p_mid = 0.5 * (p_0 + p_1);

                const std::vector<double> t_vec = {t_0, t_mid, t_1};
                const std::vector<double> p_vec = {p_0, p_mid, p_1};

                Eigen::MatrixXd points(3, 8);

                int counter = 0;
                for (int ii = 0; ii < 3; ii++) {
                    for (int jj = 0; jj < 3; jj++) {
                        if ((ii != 1) || (jj != 1)) {
                            points(0, counter) = cos(p_vec[ii]) * sin(t_vec[jj]);
                            points(1, counter) = sin(p_vec[ii]) * sin(t_vec[jj]);
                            points(2, counter) = cos(t_vec[jj]);
                            counter++;
                        }
                    }
                }

                Eigen::MatrixXd rot_1(3, 3);
                Eigen::MatrixXd rot_2(3, 3);

                rot_1(0, 0) = cos(p_mid);
                rot_1(0, 1) = sin(p_mid);
                rot_1(0, 2) = 0.0;
                rot_1(1, 0) = -sin(p_mid);
                rot_1(1, 1) = cos(p_mid);
                rot_1(1, 2) = 0.0;
                rot_1(2, 0) = 0.0;
                rot_1(2, 1) = 0.0;
                rot_1(2, 2) = 1.0;

                rot_2(0, 0) = sin(t_mid);
                rot_2(0, 1) = 0.0;
                rot_2(0, 2) = cos(t_mid);
                rot_2(1, 0) = 0.0;
                rot_2(1, 1) = 1.0;
                rot_2(1, 2) = 0.0;
                rot_2(2, 0) = cos(t_mid);
                rot_2(2, 1) = 0.0;
                rot_2(2, 2) = -sin(t_mid);

                Eigen::MatrixXd rotated_points = rot_2 * rot_1 * points;

                double t_min = std::numeric_limits<double>::max(), t_max = std::numeric_limits<double>::lowest();
                double p_min = std::numeric_limits<double>::max(), p_max = std::numeric_limits<double>::lowest();

                for (int ii = 0; ii < 8; ii++) {

                    const double rho = std::sqrt(rotated_points(0,ii)*rotated_points(0,ii) + rotated_points(1,ii)*rotated_points(1,ii) + rotated_points(2,ii)*rotated_points(2,ii));
                    const double phi = ((rotated_points(1,ii) > 0) ? 1 : ((rotated_points(1,ii) < 0) ? -1 : 0)) * acos(rotated_points(0,ii) / sqrt(rotated_points(0,ii)*rotated_points(0,ii) + rotated_points(1,ii)*rotated_points(1,ii)));
                    const double theta = acos(rotated_points(2,ii) / rho);

                    if (t_min > theta) t_min = theta;
                    if (t_max < theta) t_max = theta;
                    if (p_min > phi) p_min = phi;
                    if (p_max < phi) p_max = phi;

                }

                const double size_theta = t_max - t_min;
                const double size_phi = p_max - p_min;

                t_min -= 0.05 * size_theta;
                t_max += 0.05 * size_theta;
                p_min -= 0.05 * size_phi;
                p_max += 0.05 * size_phi;

                std::vector<double> vec_phi(PS_II, 0.0);
                for (long long jj = 0; jj < PS_II; jj++) {
                    if (TYPE_POINTS == 0) {
                        vec_phi[jj] = 0.5 * (p_max - p_min) * std::cos(M_PI * (2.0 * jj + 1.0) / (2.0 * PS_II)) + 0.5 * (p_max + p_min);
                    } else {
                        vec_phi[jj] = p_min + jj * (p_max - p_min) / (PS_II - 1);
                    }
                }

                std::vector<double> vec_theta(PT_II, 0.0);
                for (long long jj = 0; jj < PT_II; jj++) {
                    if (TYPE_POINTS == 0) {
                        vec_theta[jj] = 0.5 * (t_max - t_min) * std::cos(M_PI * (2.0 * jj + 1.0) / (2.0 * PT_II)) + 0.5 * (t_max + t_min);
                    } else {
                        vec_theta[jj] = t_min + jj * (t_max - t_min) / (PT_II - 1);
                    }
                }

                Eigen::MatrixXd new_points(3, PS_II*PT_II);
                    
                for (long long ii = 0; ii < PS_II; ii++) {
                    for (long long jj = 0; jj < PT_II; jj++) {
                        new_points(0, ii*PT_II+jj) = cos(vec_phi[ii]) * sin(vec_theta[jj]);
                        new_points(1, ii*PT_II+jj) = sin(vec_phi[ii]) * sin(vec_theta[jj]);
                        new_points(2, ii*PT_II+jj) = cos(vec_theta[jj]);
                    }
                }

                new_points = rot_1.inverse() * rot_2.inverse() * new_points; 

                interpolation_points_thread[nonrel_conesegment] = new_points;

            }

            for (const auto & [nonrel_conesegment, new_points] : interpolation_points_thread) {

                #pragma omp critical

                levels_[level]->interpolation_points_[nonrel_conesegment] = new_points;

            }

            }

        }

        // template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        // const double, const std::complex<double>, const std::complex<double>, std::complex<double>&)>
        // void LevelDIndexes()

    void LevelDIndexes(
        const std::function<void(
        double, double, double, double, double, double, double, double, double, std::complex<double>,
        std::complex<double>, std::complex<double>, std::complex<double>&)>& _kernel)
        {

            #pragma omp parallel
            {

            long long old_cocentered_morton_box = -1;

            long long npoints;
            long long points_begin;
            std::array<double, 3> center;

            #pragma omp for
            for (long long i = 0; i < levels_.back()->relconesegment2nonrelconesegment_.size(); i++) {

                const long long cocentered_morton_box = levels_.back()->relconesegment2cocenteredmortonboxid_[i];
                const long long nonrel_conesegment = levels_.back()->relconesegment2nonrelconesegment_[i];
                const double radius1 = levels_.back()->relconesegment2minradius_[i];
                const double radius2 = levels_.back()->relconesegment2maxradius_[i];
                
                const long long coefficients_begin = i * P_; 
                
                if (old_cocentered_morton_box != cocentered_morton_box) {

                    const std::array<long long, 2>& points_data = levels_.back()->mortonbox2discretizationpoints_2_.at(cocentered_morton_box);
                    points_begin = points_data[0];
                    npoints = points_data[1];

                    old_cocentered_morton_box = cocentered_morton_box;

                    levels_.back()->GetBoxCenter(cocentered_morton_box, center[0], center[1], center[2]);

                }
                
                Eigen::MatrixXd new_points = levels_.back()->interpolation_points_[nonrel_conesegment];

                Eigen::MatrixXcd A(2*PS_II*PT_II, npoints);               

                for (long long i1 = 0; i1 < PS_II; i1++) {
                    for (long long i2 = 0; i2 < PT_II; i2++) {

                        const double x = center[0] + radius1 * new_points(0, i1*PT_II+i2);
                        const double y = center[1] + radius1 * new_points(1, i1*PT_II+i2);
                        const double z = center[2] + radius1 * new_points(2, i1*PT_II+i2);

                        for (int jj = 0; jj < npoints; jj++) {

                            std::complex<double> kernel;   

                            const double locx = x_[points_begin + jj];
                            const double locy = y_[points_begin + jj];
                            const double locz = z_[points_begin + jj];

                            const double loc_normal_x = normal_x_[points_begin + jj];
                            const double loc_normal_y = normal_y_[points_begin + jj];
                            const double loc_normal_z = normal_z_[points_begin + jj];

                            const std::complex<double> locdensity = {1.0, 0.0};

                            _kernel(locx, locy, locz,
                                    x, y, z,
                                    loc_normal_x, loc_normal_y, loc_normal_z,
                                    coupling_parameter_, wavenumber_,
                                    locdensity,
                                    kernel);

                            A(i1 * PT_II + i2, jj) = kernel;

                        }

                    }
                }  

                for (long long i1 = 0; i1 < PS_II; i1++) {
                    for (long long i2 = 0; i2 < PT_II; i2++) {

                        const double x = center[0] + radius2 * new_points(0, i1*PT_II+i2);
                        const double y = center[1] + radius2 * new_points(1, i1*PT_II+i2);
                        const double z = center[2] + radius2 * new_points(2, i1*PT_II+i2);

                        for (int jj = 0; jj < npoints; jj++) {

                            std::complex<double> kernel;   

                            const double locx = x_[points_begin + jj];
                            const double locy = y_[points_begin + jj];
                            const double locz = z_[points_begin + jj];

                            const double loc_normal_x = normal_x_[points_begin + jj];
                            const double loc_normal_y = normal_y_[points_begin + jj];
                            const double loc_normal_z = normal_z_[points_begin + jj];

                            const std::complex<double> locdensity = {1.0, 0.0};

                            _kernel(locx, locy, locz,
                                    x, y, z,
                                    loc_normal_x, loc_normal_y, loc_normal_z,
                                    coupling_parameter_, wavenumber_,
                                    locdensity,
                                    kernel);

                            A(PS_II*PT_II + i1 * PT_II + i2, jj) = kernel;

                        }

                    }
                }  

                Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> qr = A.colPivHouseholderQr();

                long long num_idxs;

                if (TOTAL_COLS != -1) {

                    if (P_ < npoints) {
                        num_idxs = P_;
                    } else {
                        num_idxs = npoints;
                    }

                } else {

                    Eigen::VectorXcd diag = qr.matrixQR().diagonal();
                    bool condition = true;
                    long long idx = 0;
                    while (condition) {
                        const double norm = std::sqrt(diag(idx).real()*diag(idx).real() + diag(idx).imag()*diag(idx).imag());
                        idx++;
                        condition = (norm > ACCURACY) && (idx < diag.size());
                    }
                    num_idxs = idx;

                }
                
                Eigen::VectorXi P = qr.colsPermutation().indices();  

                for (long long ii = 0; ii < num_idxs; ii++) {
                    
                    levels_.back()->indexes_relconesegment_[coefficients_begin + ii] = points_begin + P(ii);

                }    

                levels_.back()->num_elems_indexes_relconesegment_[i] = num_idxs;

            }

            }
 
        }

        void CommunicateIndexes(const int level,
                                std::unordered_map<long long, std::vector<long long>> & indexes) const 
        {

            indexes.clear();

            for (int rank = 0; rank < comm_size_; rank++) {

                std::vector<long long> num_idxs_loc;
                std::vector<long long> num_idxs_all;

                std::vector<long long> indexes_loc;
                std::vector<long long> indexes_all;

                std::vector<int> recv_counts_indexes;
                std::vector<int> displs_indexes;

                std::vector<int> recv_counts_num_idxs;
                std::vector<int> displs_num_idxs;

                if (rank == comm_rank_) {

                    for (const auto & cs : levels_[level]->cone_segments_propagation_not_in_rank_all_) {

                        const long long idx = (cs - levels_[level]->split_points_relconesegments_[rank]) * P_;

                        num_idxs_all.push_back(levels_[level]->num_elems_indexes_relconesegment_[idx / P_]);
                        indexes_all.insert(indexes_all.end(), levels_[level]->indexes_relconesegment_.begin() + idx, levels_[level]->indexes_relconesegment_.begin() + idx + P_);

                    }

                    recv_counts_indexes = levels_[level]->recv_counts_propagation_;
                    displs_indexes = levels_[level]->displs_propagation_;

                    recv_counts_num_idxs = levels_[level]->recv_counts_propagation_;
                    displs_num_idxs = levels_[level]->displs_propagation_;

                    for (int i = 0; i < comm_size_; i++) {
                        recv_counts_num_idxs[i] /= P_;
                        displs_num_idxs[i] /= P_;
                    }

                }

                int size_loc = levels_[level]->propagation_not_in_rank_[rank].size() * P_;
                
                num_idxs_loc.resize(size_loc / P_);
                indexes_loc.resize(size_loc);
                
                MPI_Scatterv(&indexes_all[0], &recv_counts_indexes[0], &displs_indexes[0], MPI_LONG_LONG, &indexes_loc[0], size_loc, MPI_LONG_LONG, rank, mpi_comm_);
                MPI_Scatterv(&num_idxs_all[0], &recv_counts_num_idxs[0], &displs_num_idxs[0], MPI_LONG_LONG, &num_idxs_loc[0], size_loc/P_, MPI_LONG_LONG, rank, mpi_comm_);

                for (long long i = 0; i < levels_[level]->propagation_not_in_rank_[rank].size(); i++) {

                    const long long cone_segment = levels_[level]->propagation_not_in_rank_[rank][i];

                    const long long num_indx = num_idxs_loc[i];
                    const std::vector<long long> indexes_vec(indexes_loc.begin() + i * P_, indexes_loc.begin() + i * P_ + num_indx);
                    
                    indexes[cone_segment] = indexes_vec;

                }

            }

        }

        // template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        //     const double, const std::complex<double>, const std::complex<double>, std::complex<double>&)>
        // void LeveldIndexes(int level,
        //                   const std::unordered_map<long long, std::vector<long long>> & indexes_not_in_rank)
        void LeveldIndexes(
    int level,
    const std::unordered_map<long long, std::vector<long long>>& indexes_not_in_rank,
    const std::function<void(
        double, double, double, double, double, double, double, double, double, std::complex<double>,
        std::complex<double>, std::complex<double>, std::complex<double>&)>& _kernel)
        {

            #pragma omp parallel
            {

            long long old_cocentered_morton_box = -1;
            std::vector<long long> morton_children;
            std::array<double, 3> center;

            #pragma omp for
            for (long long i = 0; i < levels_[level-1]->relconesegment2nonrelconesegment_.size(); i++) {

                const long long cocentered_morton_box = levels_[level-1]->relconesegment2cocenteredmortonboxid_[i];
                const long long nonrel_conesegment = levels_[level-1]->relconesegment2nonrelconesegment_[i];
                const double radius1 = levels_[level-1]->relconesegment2minradius_[i];
                const double radius2 = levels_[level-1]->relconesegment2maxradius_[i];

                Eigen::MatrixXd new_points = levels_[level-1]->interpolation_points_[nonrel_conesegment];

                const long long new_coeffs_begin_id = i * P_; 
                
                if (old_cocentered_morton_box != cocentered_morton_box) {

                    morton_children = GetMortonRelChildren(level-1, cocentered_morton_box);

                    old_cocentered_morton_box = cocentered_morton_box;

                    levels_[level-1]->GetBoxCenter(cocentered_morton_box, center[0], center[1], center[2]);

                }

                std::unordered_set<long long> relevant_points_set;

                for (int childiter = 0; childiter < morton_children.size(); childiter++) {

                    const long long morton_child = morton_children[childiter];

                    const std::array<long long, 2> data_child = levels_[level]->mortonbox2discretizationpoints_2_[morton_child];

                    if (USE_ADAPTIVITY && (data_child[1] <= MAX_ELEMS_LEAF)) {

                        for (long long ii = 0; ii < data_child[1]; ii++) {

                            relevant_points_set.insert(data_child[0] + ii);

                        }

                    } else {

                        std::unordered_set<long long> rel_conesegments;

                        for (int interpiter = 0; interpiter < PT_II*PS_II; interpiter++) {

                            const double x = center[0] + radius1 * new_points(0, interpiter);
                            const double y = center[1] + radius1 * new_points(1, interpiter);
                            const double z = center[2] + radius1 * new_points(2, interpiter);

                            const long long locnonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_child, x, y, z);

                            const long long locrel_conesegment = levels_[level]->mortonboxnonrelcone2relcone_2_.at(morton_child).at(locnonrel_conesegment);
                            
                            rel_conesegments.insert(locrel_conesegment);

                        }

                        for (int interpiter = 0; interpiter < PT_II*PS_II; interpiter++) {

                            const double x = center[0] + radius2 * new_points(0, interpiter);
                            const double y = center[1] + radius2 * new_points(1, interpiter);
                            const double z = center[2] + radius2 * new_points(2, interpiter);

                            const long long locnonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_child, x, y, z);

                            const long long locrel_conesegment = levels_[level]->mortonboxnonrelcone2relcone_2_.at(morton_child).at(locnonrel_conesegment);
                            
                            rel_conesegments.insert(locrel_conesegment);

                        }

                        for (long long rel_cone : rel_conesegments) {

                            if ((levels_[level]->split_points_relconesegments_[comm_rank_] <= rel_cone) && (levels_[level]->split_points_relconesegments_[comm_rank_+1] > rel_cone)) {

                                const long long coeffs_begin_id = (rel_cone - levels_[level]->split_points_relconesegments_[comm_rank_]) * P_; 
                                const long long num_idx = levels_[level]->num_elems_indexes_relconesegment_[coeffs_begin_id / P_];
                                relevant_points_set.insert(levels_[level]->indexes_relconesegment_.begin() + coeffs_begin_id, levels_[level]->indexes_relconesegment_.begin() + coeffs_begin_id + num_idx);
                                
                            } else {
                                
                                const std::vector<long long> indexes = indexes_not_in_rank.at(rel_cone);
                                relevant_points_set.insert(indexes.begin(), indexes.end());
                                
                            }

                        }

                    }
                
                }

                std::vector<long long> relevant_points(relevant_points_set.begin(), relevant_points_set.end());
                relevant_points_set.clear();
                long long npoints = relevant_points.size();

                Eigen::MatrixXcd A(2*PS_II*PT_II, npoints);
                
                for (long long i1 = 0; i1 < PS_II; i1++) {
                    for (long long i2 = 0; i2 < PT_II; i2++) {

                        const double x = center[0] + radius1 * new_points(0, i1*PT_II+i2);
                        const double y = center[1] + radius1 * new_points(1, i1*PT_II+i2);
                        const double z = center[2] + radius1 * new_points(2, i1*PT_II+i2);

                        for (long long jj = 0; jj < npoints; jj++) {

                            std::complex<double> kernel;   

                            const double locx = x_[relevant_points[jj]];
                            const double locy = y_[relevant_points[jj]];
                            const double locz = z_[relevant_points[jj]];

                            const double loc_normal_x = normal_x_[relevant_points[jj]];
                            const double loc_normal_y = normal_y_[relevant_points[jj]];
                            const double loc_normal_z = normal_z_[relevant_points[jj]];

                            const std::complex<double> locdensity = {1.0, 0.0};

                            _kernel(locx, locy, locz,
                                    x, y, z,
                                    loc_normal_x, loc_normal_y, loc_normal_z,
                                    coupling_parameter_, wavenumber_,
                                    locdensity,
                                    kernel);

                            A(i1 * PT_II + i2, jj) = kernel;

                        }

                    }
                }  

                for (long long i1 = 0; i1 < PS_II; i1++) {
                    for (long long i2 = 0; i2 < PT_II; i2++) {

                        const double x = center[0] + radius2 * new_points(0, i1*PT_II+i2);
                        const double y = center[1] + radius2 * new_points(1, i1*PT_II+i2);
                        const double z = center[2] + radius2 * new_points(2, i1*PT_II+i2);

                        for (long long jj = 0; jj < npoints; jj++) {

                            std::complex<double> kernel;   

                            const double locx = x_[relevant_points[jj]];
                            const double locy = y_[relevant_points[jj]];
                            const double locz = z_[relevant_points[jj]];

                            const double loc_normal_x = normal_x_[relevant_points[jj]];
                            const double loc_normal_y = normal_y_[relevant_points[jj]];
                            const double loc_normal_z = normal_z_[relevant_points[jj]];

                            const std::complex<double> locdensity = {1.0, 0.0};

                            _kernel(locx, locy, locz,
                                    x, y, z,
                                    loc_normal_x, loc_normal_y, loc_normal_z,
                                    coupling_parameter_, wavenumber_,
                                    locdensity,
                                    kernel);

                            A(PS_II*PT_II + i1 * PT_II + i2, jj) = kernel;

                        }

                    }
                }  

                Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> qr = A.colPivHouseholderQr();

                long long num_idxs;

                if (TOTAL_COLS != -1) {

                    if (P_ < A.cols()) {
                        num_idxs = P_;
                    } else {
                        num_idxs = A.cols();
                    }

                } else {

                    Eigen::VectorXcd diag = qr.matrixQR().diagonal();
                    bool condition = true;
                    long long idx = 0;
                    while (condition) {
                        const double norm = std::sqrt(diag(idx).real()*diag(idx).real() + diag(idx).imag()*diag(idx).imag());
                        idx++;
                        condition = (norm > ACCURACY) && (idx < diag.size());
                    }
                    num_idxs = idx;

                }

                Eigen::VectorXi P = qr.colsPermutation().indices();  

                for (long long ii = 0; ii < num_idxs; ii++) {
                    
                    levels_[level-1]->indexes_relconesegment_[new_coeffs_begin_id + ii] = relevant_points[P(ii)];

                }    

                levels_[level-1]->num_elems_indexes_relconesegment_[i] = num_idxs;

            }

            }
 
        }

        void ZeroSolution() 
        {

            #pragma omp parallel for
            for (long long i = 0; i < N_loc_orig_; i++) {

                solution_[i] = {0.0, 0.0};

            }

        }

        
        // template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        //                        const double, const std::complex<double>, const std::complex<double>, 
        //                        std::complex<double>&)>
        // void SingularInteractions(const std::vector<std::complex<double>>& density)

        void SingularInteractions(
    const std::vector<std::complex<double>>& density,
    const std::function<void(
        double, double, double, double, double, double, double, double, double, std::complex<double>,
        std::complex<double>, std::complex<double>, std::complex<double>&)>& _kernel)

        {

            #pragma omp parallel
            {
            
            long long old_mortonbox = -1;
            std::vector<long long> neighbours;
            long long npoints = 0;
            long long points_begin = 0;

            #pragma omp for
            for (long long i = 0; i < N_loc_orig_; i++) {

                const long long point = sorting_[point_low_ + i];

                const double x = x_[point_low_ + i];
                const double y = y_[point_low_ + i];
                const double z = z_[point_low_ + i];

                const long long mortonbox = levels_.back()->Point2Morton(x, y, z);
                
                std::complex<double> result = solution_[i];
                std::complex<double> tmp_result;     

                if (mortonbox != old_mortonbox) {

                    neighbours = levels_.back()->GetNeighbours(mortonbox);

                    if (neighbours.size() > 27)
                        throw std::logic_error("Cannot have more than 27 neighbours.");

                    old_mortonbox = mortonbox;

                }

                for (int neighbouriter = 0; neighbouriter < neighbours.size(); neighbouriter++) {

                    const long long neighbourmorton = neighbours[neighbouriter];

                    const std::array<long long, 2> points_data = levels_.back()->mortonbox2discretizationpoints_2_[neighbourmorton];

                    points_begin = points_data[0];
                    npoints = points_data[1];

                    for (long long sourceiter = 0; sourceiter < npoints; sourceiter++) {

                        const long long num_patch_source = sorting_[points_begin + sourceiter] / (n_pts_per_patch_[0] * n_pts_per_patch_[1]);
                        bool found = precomputations_data_[point].count(num_patch_source) != 0;
                        //bool found = (sourceiter+points_begin == point_low_+i);

                        if (found) {
                            continue;
                        }

                        const double locx = x_[points_begin + sourceiter];
                        const double locy = y_[points_begin + sourceiter];
                        const double locz = z_[points_begin + sourceiter];

                        const double loc_normal_x = normal_x_[points_begin + sourceiter];
                        const double loc_normal_y = normal_y_[points_begin + sourceiter];
                        const double loc_normal_z = normal_z_[points_begin + sourceiter];

                        const std::complex<double> locdensity = density[points_begin + sourceiter];

                        _kernel(locx, locy, locz, x, y, z, loc_normal_x, loc_normal_y, loc_normal_z, coupling_parameter_, wavenumber_, locdensity, tmp_result);
                        
                        result += tmp_result;

                    }    

                }

                solution_[i] = result;

            }
            
            }

        }

        // template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        //     const double, const std::complex<double>, const std::complex<double>, std::complex<double>&),
        //     void _factorization(const double, const std::complex<double>, std::complex<double>&)>
        // void LevelDEvaluations(const std::vector<std::complex<double>>& density,
        //                        std::vector<std::complex<double>>& conesegments_current,
        //                        std::vector<std::complex<double>>& conesegments_prev)
        
void LevelDEvaluations(
    const std::vector<std::complex<double>>& density,
    std::vector<std::complex<double>>& conesegments_current,
    std::vector<std::complex<double>>& conesegments_prev,
    const std::function<void(
        double, double, double, double, double, double, double, double, double, std::complex<double>,
        std::complex<double>, std::complex<double>, std::complex<double>&)>& _kernel,
    const std::function<void(
        double, std::complex<double>, std::complex<double>&)>& _factorization)
        {

            #pragma omp parallel
            {

            long long old_cocentered_morton_box = -1;
            std::complex<double> locdensity;
            double locx;
            double locy;
            double locz;
            double loc_normal_x;
            double loc_normal_y;
            double loc_normal_z;

            long long npoints;
            long long points_begin;

            #pragma omp for
            for (long long i = 0; i < levels_.back()->relconesegment2nonrelconesegment_.size(); i++) {

                const long long cocentered_morton_box = levels_.back()->relconesegment2cocenteredmortonboxid_[i];
                const long long nonrel_conesegment = levels_.back()->relconesegment2nonrelconesegment_[i];
                double x, y, z;
                const long long coefficients_begin = i * IPSCHEME_.GetNInterpPoints(); 
                
                if (old_cocentered_morton_box != cocentered_morton_box) {

                    const std::array<long long, 2>& points_data = levels_.back()->mortonbox2discretizationpoints_2_.at(cocentered_morton_box);
                    points_begin = points_data[0];
                    npoints = points_data[1];

                    old_cocentered_morton_box = cocentered_morton_box;

                }

                for (int j = 0; j < IPSCHEME_.GetNInterpPoints(); j++) {

                    IPSCHEME_.GetInterpolationPoint(j, x, y, z);

                    const double radius = levels_.back()->Cheb2Radius(nonrel_conesegment, x);

                    levels_.back()->Cheb2Cart(cocentered_morton_box, nonrel_conesegment, x, y, z);
                    
                    std::complex<double> result = {0.0, 0.0};
                    std::complex<double> kernel;

                    for (long long k = 0; k < npoints; k++) {

                        locx = x_[points_begin + k];
                        locy = y_[points_begin + k];
                        locz = z_[points_begin + k];

                        loc_normal_x = normal_x_[points_begin + k];
                        loc_normal_y = normal_y_[points_begin + k];
                        loc_normal_z = normal_z_[points_begin + k];

                        locdensity = density[points_begin + k];

                        _kernel(locx, locy, locz,
                                x, y, z,
                                loc_normal_x, loc_normal_y, loc_normal_z,
                                coupling_parameter_, wavenumber_,
                                locdensity,
                                kernel);

                        result += kernel;

                    }

                    std::complex<double> fac;
                    _factorization(radius, wavenumber_, fac);
                    const double dabsfac = 1.0/(fac.real()*fac.real() + fac.imag()*fac.imag());                    

                    std::complex<double> new_fac = {1.0, 0.0};

                    const double tmp_fac_real = fac.real() * new_fac.real() + fac.imag() * new_fac.imag();
                    const double tmp_fac_imag = fac.real() * new_fac.imag() - fac.imag() * new_fac.real();

                    fac = {tmp_fac_real, tmp_fac_imag};

                    const long long coefficients_pos = coefficients_begin + j;
                    const double real_part = (result.real() * fac.real() - result.imag() * fac.imag()) * dabsfac;
                    const double imag_part = (result.imag() * fac.real() + result.real() * fac.imag()) * dabsfac;

                    conesegments_prev[coefficients_pos] = {real_part, imag_part};

                }

                IPSCHEME_.GenerateInterpolant(&conesegments_prev[coefficients_begin]);                 

            }

            #pragma omp for
            for (long long i = 0; i < conesegments_current.size(); i++) {

                conesegments_current[i] = conesegments_prev[i];

            }

            }
            
        }        

        // template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        //     const double, const std::complex<double>, const std::complex<double>, std::complex<double>&)>
        // void LevelDEvaluations(const std::vector<std::complex<double>>& density,
        //                        std::vector<std::complex<double>>& conesegments_current,
        //                        std::vector<std::complex<double>>& conesegments_prev)
        void LevelDEvaluations(
    const std::vector<std::complex<double>>& density,
    std::vector<std::complex<double>>& conesegments_current,
    std::vector<std::complex<double>>& conesegments_prev,
    const std::function<void(
        double, double, double, double, double, double, double, double, double, std::complex<double>,
        std::complex<double>, std::complex<double>, std::complex<double>&)>& _kernel)
        {

            #pragma omp parallel
            {

            long long old_cocentered_morton_box = -1;

            long long npoints;
            long long points_begin;
            std::array<double, 3> center;

            #pragma omp for
            for (long long i = 0; i < levels_.back()->relconesegment2nonrelconesegment_.size(); i++) {

                const long long cocentered_morton_box = levels_.back()->relconesegment2cocenteredmortonboxid_[i];
                const long long nonrel_conesegment = levels_.back()->relconesegment2nonrelconesegment_[i];
                const double radius1 = levels_.back()->relconesegment2minradius_[i];
                const double radius2 = levels_.back()->relconesegment2maxradius_[i];

                double x, y, z;
                const long long coefficients_begin = i * P_;
                
                if (old_cocentered_morton_box != cocentered_morton_box) {

                    const std::array<long long, 2>& points_data = levels_.back()->mortonbox2discretizationpoints_2_.at(cocentered_morton_box);
                    points_begin = points_data[0];
                    npoints = points_data[1];

                    old_cocentered_morton_box = cocentered_morton_box;

                    levels_.back()->GetBoxCenter(cocentered_morton_box, center[0], center[1], center[2]);

                }
                
                Eigen::MatrixXd new_points = levels_.back()->interpolation_points_[nonrel_conesegment];

                Eigen::VectorXcd b(2*PS_II*PT_II);
                b.setZero();               

                for (long long i1 = 0; i1 < PS_II; i1++) {
                    for (long long i2 = 0; i2 < PT_II; i2++) {

                        const double x = center[0] + radius1 * new_points(0, i1*PT_II+i2);
                        const double y = center[1] + radius1 * new_points(1, i1*PT_II+i2);
                        const double z = center[2] + radius1 * new_points(2, i1*PT_II+i2);

                        for (long long jj = 0; jj < npoints; jj++) {

                            std::complex<double> kernel;   

                            const double locx = x_[points_begin + jj];
                            const double locy = y_[points_begin + jj];
                            const double locz = z_[points_begin + jj];

                            const double loc_normal_x = normal_x_[points_begin + jj];
                            const double loc_normal_y = normal_y_[points_begin + jj];
                            const double loc_normal_z = normal_z_[points_begin + jj];

                            const std::complex<double> locdensity = density[points_begin + jj];

                            _kernel(locx, locy, locz,
                                    x, y, z,
                                    loc_normal_x, loc_normal_y, loc_normal_z,
                                    coupling_parameter_, wavenumber_,
                                    locdensity,
                                    kernel);

                            b[i1 * PT_II + i2] += kernel;

                        }

                    }
                } 

                for (long long i1 = 0; i1 < PS_II; i1++) {
                    for (long long i2 = 0; i2 < PT_II; i2++) {

                        const double x = center[0] + radius2 * new_points(0, i1*PT_II+i2);
                        const double y = center[1] + radius2 * new_points(1, i1*PT_II+i2);
                        const double z = center[2] + radius2 * new_points(2, i1*PT_II+i2);

                        for (long long jj = 0; jj < npoints; jj++) {

                            std::complex<double> kernel;   

                            const double locx = x_[points_begin + jj];
                            const double locy = y_[points_begin + jj];
                            const double locz = z_[points_begin + jj];

                            const double loc_normal_x = normal_x_[points_begin + jj];
                            const double loc_normal_y = normal_y_[points_begin + jj];
                            const double loc_normal_z = normal_z_[points_begin + jj];

                            const std::complex<double> locdensity = density[points_begin + jj];

                            _kernel(locx, locy, locz,
                                    x, y, z,
                                    loc_normal_x, loc_normal_y, loc_normal_z,
                                    coupling_parameter_, wavenumber_,
                                    locdensity,
                                    kernel);

                            b[PS_II*PT_II + i1 * PT_II + i2] += kernel;

                        }

                    }
                } 

                const int num_idxs = levels_.back()->num_elems_indexes_relconesegment_[i];

                Eigen::MatrixXcd A(2*PS_II*PT_II, num_idxs);

                for (long long i1 = 0; i1 < PS_II; i1++) {
                    for (long long i2 = 0; i2 < PT_II; i2++) {

                        const double x = center[0] + radius1 * new_points(0, i1*PT_II+i2);
                        const double y = center[1] + radius1 * new_points(1, i1*PT_II+i2);
                        const double z = center[2] + radius1 * new_points(2, i1*PT_II+i2);

                        for (long long jj = 0; jj < num_idxs; jj++) {

                            const long long idx = levels_.back()->indexes_relconesegment_[coefficients_begin + jj];

                            std::complex<double> kernel;   

                            const double locx = x_[idx];
                            const double locy = y_[idx];
                            const double locz = z_[idx];

                            const double loc_normal_x = normal_x_[idx];
                            const double loc_normal_y = normal_y_[idx];
                            const double loc_normal_z = normal_z_[idx];

                            const std::complex<double> locdensity = {1.0, 0.0};

                            _kernel(locx, locy, locz,
                                    x, y, z,
                                    loc_normal_x, loc_normal_y, loc_normal_z,
                                    coupling_parameter_, wavenumber_,
                                    locdensity,
                                    kernel);

                            A(i1 * PT_II + i2, jj) = kernel;

                        }

                    }
                } 

                for (long long i1 = 0; i1 < PS_II; i1++) {
                    for (long long i2 = 0; i2 < PT_II; i2++) {

                        const double x = center[0] + radius2 * new_points(0, i1*PT_II+i2);
                        const double y = center[1] + radius2 * new_points(1, i1*PT_II+i2);
                        const double z = center[2] + radius2 * new_points(2, i1*PT_II+i2);

                        for (long long jj = 0; jj < num_idxs; jj++) {

                            const long long idx = levels_.back()->indexes_relconesegment_[coefficients_begin + jj];

                            std::complex<double> kernel;   

                            const double locx = x_[idx];
                            const double locy = y_[idx];
                            const double locz = z_[idx];

                            const double loc_normal_x = normal_x_[idx];
                            const double loc_normal_y = normal_y_[idx];
                            const double loc_normal_z = normal_z_[idx];

                            const std::complex<double> locdensity = {1.0, 0.0};

                            _kernel(locx, locy, locz,
                                    x, y, z,
                                    loc_normal_x, loc_normal_y, loc_normal_z,
                                    coupling_parameter_, wavenumber_,
                                    locdensity,
                                    kernel);

                            A(PS_II*PT_II + i1 * PT_II + i2, jj) = kernel;

                        }

                    }
                } 

                Eigen::VectorXcd coeffs = A.completeOrthogonalDecomposition().solve(b);
                
                for (long long ii = 0; ii < num_idxs; ii++) {

                    conesegments_prev[coefficients_begin + ii] = coeffs[ii];

                }

            }

            #pragma omp for
            for (long long i = 0; i < conesegments_current.size(); i++) {

                conesegments_current[i] = conesegments_prev[i];

            }

            }
 
        }

        void CommunicateInterpolation(const int level, std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
                                      const std::vector<std::complex<double>>& conesegments_prev) const 
        {

            coeffs.clear(); 

            for (int rank = 0; rank < comm_size_; rank++) {

                std::vector<std::complex<double>> coeffs_cone_segments_all;
                std::vector<std::complex<double>> coeffs_cone_segments_loc;

                std::vector<int> recv_counts;
                std::vector<int> displs;

                if (rank == comm_rank_) {

                    for (const auto & cs : levels_[level]->cone_segments_interpolation_not_in_rank_all_) {

                        const long long idx = (cs - levels_[level]->split_points_relconesegments_[rank]) * P_;

                        coeffs_cone_segments_all.insert(coeffs_cone_segments_all.end(), conesegments_prev.begin() + idx, conesegments_prev.begin() + idx + P_);

                    }

                    recv_counts = levels_[level]->recv_counts_interpolation_;
                    displs = levels_[level]->displs_interpolation_;

                }

                int size_loc = levels_[level]->interpolation_not_in_rank_[rank].size() * P_;
                coeffs_cone_segments_loc.resize(size_loc);

                MPI_Scatterv(&coeffs_cone_segments_all[0], &recv_counts[0], &displs[0], MPI_DOUBLE_COMPLEX, &coeffs_cone_segments_loc[0], size_loc, MPI_DOUBLE_COMPLEX, rank, mpi_comm_);

                for (long long i = 0; i < levels_[level]->interpolation_not_in_rank_[rank].size(); i++) {

                    const long long cone_segment = levels_[level]->interpolation_not_in_rank_[rank][i];

                    const std::vector<std::complex<double>> coefficients(coeffs_cone_segments_loc.begin() + i * P_, coeffs_cone_segments_loc.begin() + (i+1) * P_);

                    coeffs[cone_segment] = coefficients;

                }

            }

        }

        // template<void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        // const double, const std::complex<double>, const std::complex<double>, std::complex<double>&),
        // void _factorization(const double, const std::complex<double>, std::complex<double>&)>
        // void Interpolation(const int level, const std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
        //                    const std::vector<std::complex<double>>& density, const std::vector<std::complex<double>>& conesegments_prev)
        void Interpolation(
            int level,
            const std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
            const std::vector<std::complex<double>>& density,
            const std::vector<std::complex<double>>& conesegments_prev,
            const std::function<void(
                double, double, double, double, double, double, double, double, double, std::complex<double>,
                std::complex<double>, std::complex<double>, std::complex<double>&)>& _kernel,
            const std::function<void(
                double, std::complex<double>, std::complex<double>&)>& _factorization)
                    
        {

            #pragma omp parallel
            {
                
            long long old_morton_box = -1;
            std::vector<long long> morton_cousinboxes;            

            #pragma omp for
            for (long long pointiter = 0; pointiter < N_loc_orig_; pointiter++) {

                const double x = x_[pointiter + point_low_];
                const double y = y_[pointiter + point_low_];
                const double z = z_[pointiter + point_low_];

                const long long morton_box = levels_[level]->Point2Morton(x, y, z);

                std::complex<double> value = solution_[pointiter];  

                if (morton_box != old_morton_box) {

                    morton_cousinboxes = levels_[level]->GetCousins(morton_box);
                    old_morton_box = morton_box;

                }

                for (int cousiniter = 0; cousiniter < morton_cousinboxes.size(); cousiniter++) {

                    double locx = x;
                    double locy = y;
                    double locz = z;

                    const long long morton_cousinbox = morton_cousinboxes[cousiniter];
                    
                    const std::array<long long, 2> points_data = levels_[level]->mortonbox2discretizationpoints_2_[morton_cousinbox];
                    long long points_begin = points_data[0];
                    long long npoints = points_data[1];

                    if (USE_ADAPTIVITY && npoints <= MAX_ELEMS_LEAF) {

                        for (long long sourceiter = 0; sourceiter < npoints; sourceiter++) {

                            const double locx = x_[points_begin + sourceiter];
                            const double locy = y_[points_begin + sourceiter];
                            const double locz = z_[points_begin + sourceiter];

                            const double loc_normal_x = normal_x_[points_begin + sourceiter];
                            const double loc_normal_y = normal_y_[points_begin + sourceiter];
                            const double loc_normal_z = normal_z_[points_begin + sourceiter];

                            const std::complex<double> locdensity = density[points_begin + sourceiter];

                            std::complex<double> tmp_result;
                            _kernel(locx, locy, locz, x, y, z, loc_normal_x, loc_normal_y, loc_normal_z, coupling_parameter_, wavenumber_, locdensity, tmp_result);
                            
                            value += tmp_result;

                        }    

                    } else {

                        long long nonrelconesegment;
                        levels_[level]->Cart2Cheb(morton_cousinbox, locx, locy, locz, nonrelconesegment);
                        const long long relconesegment = levels_[level]->mortonboxnonrelcone2relcone_2_.at(morton_cousinbox).at(nonrelconesegment);
                        const double radius = levels_[level]->Cheb2Radius(nonrelconesegment, locx);
                        std::complex<double> fac;
                        _factorization(radius, wavenumber_, fac); 

                        std::complex<double> new_fac = {1.0, 0.0};   

                        const double tmp_fac_real = (fac.real() * new_fac.real() + fac.imag() * new_fac.imag()) / (new_fac.real()*new_fac.real() + new_fac.imag()*new_fac.imag());
                        const double tmp_fac_imag = (- fac.real() * new_fac.imag() + fac.imag() * new_fac.real()) / (new_fac.real()*new_fac.real() + new_fac.imag()*new_fac.imag());

                        fac = {tmp_fac_real, tmp_fac_imag};

                        std::vector<std::complex<double>> vals_begin;

                        if ((relconesegment >= levels_[level]->split_points_relconesegments_[comm_rank_]) && (relconesegment < levels_[level]->split_points_relconesegments_[comm_rank_+1])) {

                            const long long coeffs_begin_id = (relconesegment - levels_[level]->split_points_relconesegments_[comm_rank_]) * P_; 
                            vals_begin.insert(vals_begin.begin(), conesegments_prev.begin() + coeffs_begin_id, conesegments_prev.begin() + coeffs_begin_id + P_);

                        } else {

                            vals_begin = coeffs.at(relconesegment);

                        }            

                        const std::complex<double> tmp = IPSCHEME_.Interpolate(locx, locy, locz, &vals_begin[0]);
                        
                        const double value_real = tmp.real() * fac.real() - tmp.imag() * fac.imag();
                        const double value_imag = tmp.real() * fac.imag() + tmp.imag() * fac.real();

                        value += std::complex<double>{value_real, value_imag};

                    }

                }                

                solution_[pointiter] = value;

            }

            }

        } 

        void CommunicateInterpolation(const int level, std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
                                      std::unordered_map<long long, long long> & num_idxs, std::unordered_map<long long, std::vector<long long>> & indexes,
                                      const std::vector<std::complex<double>>& conesegments_prev) const 
        {

            coeffs.clear();
            num_idxs.clear();
            indexes.clear();

            for (int rank = 0; rank < comm_size_; rank++) {

                std::vector<std::complex<double>> coeffs_cone_segments_all;
                std::vector<std::complex<double>> coeffs_cone_segments_loc;

                std::vector<long long> num_idxs_loc;
                std::vector<long long> num_idxs_all;

                std::vector<long long> indexes_loc;
                std::vector<long long> indexes_all;

                std::vector<int> recv_counts_coeffs;
                std::vector<int> displs_coeffs;

                std::vector<int> recv_counts_num_idxs;
                std::vector<int> displs_num_idxs;

                if (rank == comm_rank_) {

                    for (const auto & cs : levels_[level]->cone_segments_interpolation_not_in_rank_all_) {

                        const long long idx = (cs - levels_[level]->split_points_relconesegments_[rank]) * P_;

                        coeffs_cone_segments_all.insert(coeffs_cone_segments_all.end(), conesegments_prev.begin() + idx, conesegments_prev.begin() + idx + P_);
                        num_idxs_all.push_back(levels_[level]->num_elems_indexes_relconesegment_[idx / P_]);
                        indexes_all.insert(indexes_all.end(), levels_[level]->indexes_relconesegment_.begin() + idx, levels_[level]->indexes_relconesegment_.begin() + idx + P_);

                    }

                    recv_counts_coeffs = levels_[level]->recv_counts_interpolation_;
                    displs_coeffs= levels_[level]->displs_interpolation_;

                    recv_counts_num_idxs = levels_[level]->recv_counts_interpolation_;
                    displs_num_idxs = levels_[level]->displs_interpolation_;

                    for (int i = 0; i < comm_size_; i++) {
                        recv_counts_num_idxs[i] /= P_;
                        displs_num_idxs[i] /= P_;
                    }

                }

                int size_loc = levels_[level]->interpolation_not_in_rank_[rank].size() * P_;

                coeffs_cone_segments_loc.resize(size_loc);
                num_idxs_loc.resize(size_loc / P_);
                indexes_loc.resize(size_loc);

                MPI_Scatterv(&coeffs_cone_segments_all[0], &recv_counts_coeffs[0], &displs_coeffs[0], MPI_DOUBLE_COMPLEX, &coeffs_cone_segments_loc[0], size_loc, MPI_DOUBLE_COMPLEX, rank, mpi_comm_);
                MPI_Scatterv(&indexes_all[0], &recv_counts_coeffs[0], &displs_coeffs[0], MPI_LONG_LONG, &indexes_loc[0], size_loc, MPI_LONG_LONG, rank, mpi_comm_);
                MPI_Scatterv(&num_idxs_all[0], &recv_counts_num_idxs[0], &displs_num_idxs[0], MPI_LONG_LONG, &num_idxs_loc[0], size_loc/P_, MPI_LONG_LONG, rank, mpi_comm_);

                for (long long i = 0; i < levels_[level]->interpolation_not_in_rank_[rank].size(); i++) {

                    const long long cone_segment = levels_[level]->interpolation_not_in_rank_[rank][i];

                    const std::vector<std::complex<double>> coefficients(coeffs_cone_segments_loc.begin() + i * P_, coeffs_cone_segments_loc.begin() + (i+1) * P_);
                    const std::vector<long long> indexes_vec(indexes_loc.begin() + i * P_, indexes_loc.begin() + (i+1) * P_);
                    const long long num_indx = num_idxs_loc[i];

                    coeffs[cone_segment] = coefficients;
                    num_idxs[cone_segment] = num_indx;
                    indexes[cone_segment] = indexes_vec;

                }

            }

        }

        // template<void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        // const double, const std::complex<double>, const std::complex<double>, std::complex<double>&)>
        // void Interpolation(const int level, const std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
        //                    const std::unordered_map<long long, long long>& num_indxs, const std::unordered_map<long long, std::vector<long long>>& indexes,
        //                    const std::vector<std::complex<double>>& density,
        //                    const std::vector<std::complex<double>>& conesegments_prev)
        
        // template<void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        // const double, const std::complex<double>, const std::complex<double>, std::complex<double>&)>
        void Interpolation(const int level, const std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
                           const std::unordered_map<long long, long long>& num_indxs, const std::unordered_map<long long, std::vector<long long>>& indexes,
                           const std::vector<std::complex<double>>& density,
                           const std::vector<std::complex<double>>& conesegments_prev, 
                         const std::function<void(
                        double, double, double, double, double, double, double, double, double, std::complex<double>,
                        std::complex<double>, std::complex<double>, std::complex<double>&)>& _kernel)  
        {

            #pragma omp parallel
            {
                
            long long old_morton_box = -1;
            std::vector<long long> morton_cousinboxes;            

            #pragma omp for
            for (long long pointiter = 0; pointiter < N_loc_orig_; pointiter++) {

                const double x = x_[pointiter + point_low_];
                const double y = y_[pointiter + point_low_];
                const double z = z_[pointiter + point_low_];

                const long long morton_box = levels_[level]->Point2Morton(x, y, z);

                std::complex<double> value = solution_[pointiter];

                if (morton_box != old_morton_box) {

                    morton_cousinboxes = levels_[level]->GetCousins(morton_box);
                    old_morton_box = morton_box;

                }

                for (int cousiniter = 0; cousiniter < morton_cousinboxes.size(); cousiniter++) {

                    double locx = x;
                    double locy = y;
                    double locz = z;

                    const long long morton_cousinbox = morton_cousinboxes[cousiniter];
                    
                    const std::array<long long, 2> points_data = levels_[level]->mortonbox2discretizationpoints_2_[morton_cousinbox];
                    long long points_begin = points_data[0];
                    long long npoints = points_data[1];

                    if (USE_ADAPTIVITY && npoints <= MAX_ELEMS_LEAF) {

                        for (long long sourceiter = 0; sourceiter < npoints; sourceiter++) {

                            const double locx = x_[points_begin + sourceiter];
                            const double locy = y_[points_begin + sourceiter];
                            const double locz = z_[points_begin + sourceiter];

                            const double loc_normal_x = normal_x_[points_begin + sourceiter];
                            const double loc_normal_y = normal_y_[points_begin + sourceiter];
                            const double loc_normal_z = normal_z_[points_begin + sourceiter];

                            const std::complex<double> locdensity = density[points_begin + sourceiter];

                            std::complex<double> tmp_result;
                            _kernel(locx, locy, locz, x, y, z, loc_normal_x, loc_normal_y, loc_normal_z, coupling_parameter_, wavenumber_, locdensity, tmp_result);
                            
                            value += tmp_result;

                        }    

                    } else {

                        long long nonrelconesegment;
                        levels_[level]->Cart2Cheb(morton_cousinbox, locx, locy, locz, nonrelconesegment);
                        const long long relconesegment = levels_[level]->mortonboxnonrelcone2relcone_2_.at(morton_cousinbox).at(nonrelconesegment);
                        
                        std::vector<std::complex<double>> vals_begin;
                        long long idx_num;
                        std::vector<long long> index;

                        if ((relconesegment >= levels_[level]->split_points_relconesegments_[comm_rank_]) && (relconesegment < levels_[level]->split_points_relconesegments_[comm_rank_+1])) {

                            const long long coeffs_begin_id = (relconesegment - levels_[level]->split_points_relconesegments_[comm_rank_]) * P_; 
                            vals_begin.insert(vals_begin.begin(), conesegments_prev.begin() + coeffs_begin_id, conesegments_prev.begin() + coeffs_begin_id + P_);
                            idx_num = levels_[level]->num_elems_indexes_relconesegment_[coeffs_begin_id / P_];
                            index.insert(index.begin(), levels_[level]->indexes_relconesegment_.begin() + coeffs_begin_id, levels_[level]->indexes_relconesegment_.begin() + coeffs_begin_id + P_);

                        } else {

                            vals_begin = coeffs.at(relconesegment);
                            idx_num = num_indxs.at(relconesegment);
                            index = indexes.at(relconesegment);

                        }

                        for (long long ii = 0; ii < idx_num; ii++) {

                            const double locx = x_[index[ii]];
                            const double locy = y_[index[ii]];
                            const double locz = z_[index[ii]];

                            const double loc_normal_x = normal_x_[index[ii]];
                            const double loc_normal_y = normal_y_[index[ii]];
                            const double loc_normal_z = normal_z_[index[ii]];

                            const std::complex<double> locdensity = vals_begin[ii];

                            std::complex<double> tmp_result;
                            _kernel(locx, locy, locz, x, y, z, loc_normal_x, loc_normal_y, loc_normal_z, coupling_parameter_, wavenumber_, locdensity, tmp_result);
                            
                            value += tmp_result;

                        }

                    }

                }                

                solution_[pointiter] = value;

            }

            }

        } 

        void SwapCoefficients(std::vector<std::complex<double>>& conesegments_current,
                              std::vector<std::complex<double>>& conesegments_prev) {

            #pragma omp parallel for
            for (long long i = 0; i < conesegments_current.size(); i++) {

                std::complex<double> tmp = conesegments_current[i];
                conesegments_current[i] = conesegments_prev[i];
                conesegments_prev[i] = tmp;

            }

        }

        void CommunicatePropagation(const int level, std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
                                    const std::vector<std::complex<double>>& conesegments_current) const {

            coeffs.clear();

            for (int rank = 0; rank < comm_size_; rank++) {

                std::vector<std::complex<double>> coeffs_cone_segments_all;
                std::vector<std::complex<double>> coeffs_cone_segments_loc;

                std::vector<int> recv_counts;
                std::vector<int> displs;

                if (rank == comm_rank_) {

                    for (const auto & cs : levels_[level]->cone_segments_propagation_not_in_rank_all_) {

                        const long long idx = (cs - levels_[level]->split_points_relconesegments_[rank]) * P_;

                        coeffs_cone_segments_all.insert(coeffs_cone_segments_all.end(), conesegments_current.begin() + idx, conesegments_current.begin() + idx + P_);

                    }

                    recv_counts = levels_[level]->recv_counts_propagation_;
                    displs = levels_[level]->displs_propagation_;

                }

                int size_loc = levels_[level]->propagation_not_in_rank_[rank].size() * P_;
                coeffs_cone_segments_loc.resize(size_loc);

                MPI_Scatterv(&coeffs_cone_segments_all[0], &recv_counts[0], &displs[0], MPI_DOUBLE_COMPLEX, &coeffs_cone_segments_loc[0], size_loc, MPI_DOUBLE_COMPLEX, rank, mpi_comm_);

                for (long long i = 0; i < levels_[level]->propagation_not_in_rank_[rank].size(); i++) {

                    const long long cone_segment = levels_[level]->propagation_not_in_rank_[rank][i];

                    const std::vector<std::complex<double>> coefficients(coeffs_cone_segments_loc.begin() + i * P_, coeffs_cone_segments_loc.begin() + (i+1) * P_);

                    coeffs[cone_segment] = coefficients;

                }

            }

        }

        // template<void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        // const double, const std::complex<double>, const std::complex<double>, std::complex<double>&),
        // void _factorization(const double, const std::complex<double>, std::complex<double>&)>
        // void Propagation(const int level, const std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
        //                  const std::vector<std::complex<double>>& density,
        //                  std::vector<std::complex<double>>& conesegments_current, const std::vector<std::complex<double>>& conesegments_prev)
        void Propagation(
    int level,
    const std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
    const std::vector<std::complex<double>>& density,
    std::vector<std::complex<double>>& conesegments_current,
    const std::vector<std::complex<double>>& conesegments_prev,
    const std::function<void(
        double, double, double, double, double, double, double, double, double, std::complex<double>,
        std::complex<double>, std::complex<double>, std::complex<double>&)>& _kernel,
    const std::function<void(
        double, std::complex<double>, std::complex<double>&)>& _factorization)
        {

            #pragma omp parallel
            {

            std::vector<double> x(P_), y(P_), z(P_), radii(P_);
            std::vector<std::complex<double>> value(P_);
            double locx, locy, locz;
            long long old_cocentered_morton_box = -1;
            std::vector<long long> morton_children;
            
            #pragma omp for
            for (long long i = 0; i < levels_[level-1]->split_points_relconesegments_[comm_rank_+1]-levels_[level-1]->split_points_relconesegments_[comm_rank_]; i++) {

                const long long cocentered_morton_box = levels_[level-1]->relconesegment2cocenteredmortonboxid_[i];
                const long long nonrel_conesegment = levels_[level-1]->relconesegment2nonrelconesegment_[i];

                const long long new_coeffs_begin_id = i * P_; 
                
                if (cocentered_morton_box != old_cocentered_morton_box) {

                    morton_children = GetMortonRelChildren(level-1, cocentered_morton_box);

                    old_cocentered_morton_box = cocentered_morton_box;

                } 

                const std::array<std::unordered_map<long long, std::vector<int>>, 8>& allchildnonrelcones = levels_[level]->protobox_.at(nonrel_conesegment);

                for (int interpid = 0; interpid < IPSCHEME_.GetNInterpPoints(); interpid++) {

                    IPSCHEME_.GetInterpolationPoint(interpid, x[interpid], y[interpid], z[interpid]);
                    radii[interpid] = levels_[level-1]->Cheb2Radius(nonrel_conesegment, x[interpid]);

                    levels_[level-1]->Cheb2Cart(cocentered_morton_box, nonrel_conesegment, x[interpid], y[interpid], z[interpid]);

                    value[interpid] = {0.0, 0.0};

                }

                for (int childiter = 0; childiter < morton_children.size(); childiter++) {

                    const long long morton_child = morton_children[childiter];
                    const long long child_pos = morton_child % 8;
                    const std::unordered_map<long long, std::vector<int>>& childnonrelcones = allchildnonrelcones.at(child_pos);

                    std::array<long long, 2> data_box = levels_[level]->mortonbox2discretizationpoints_2_[morton_child];
                    const long long points_begin = data_box[0];
                    const long long npoints = data_box[1];

                    for (auto childconeiter = childnonrelcones.begin(); childconeiter != childnonrelcones.end(); childconeiter++) {

                        for (int interpiter = 0; interpiter < childconeiter->second.size(); interpiter++) {

                            const int interpid = childconeiter->second[interpiter];

                            if (USE_ADAPTIVITY && npoints <= MAX_ELEMS_LEAF) {

                                for (long long sourceiter = 0; sourceiter < npoints; sourceiter++) {

                                    const double locx = x_[points_begin + sourceiter];
                                    const double locy = y_[points_begin + sourceiter];
                                    const double locz = z_[points_begin + sourceiter];

                                    const double loc_normal_x = normal_x_[points_begin + sourceiter];
                                    const double loc_normal_y = normal_y_[points_begin + sourceiter];
                                    const double loc_normal_z = normal_z_[points_begin + sourceiter];

                                    const std::complex<double> locdensity = density[points_begin + sourceiter];

                                    std::complex<double> tmp_result;
                                    _kernel(locx, locy, locz, x[interpid], y[interpid], z[interpid], loc_normal_x, loc_normal_y, loc_normal_z, coupling_parameter_, wavenumber_, locdensity, tmp_result);
                                    
                                    value[interpid] += tmp_result;

                                } 

                            } else {
                                
                                locx = x[interpid];
                                locy = y[interpid];
                                locz = z[interpid];
                                levels_[level]->Cart2Cheb(morton_child, childconeiter->first, locx, locy, locz);

                                const long long locrelconesegment = levels_[level]->mortonboxnonrelcone2relcone_2_.at(morton_child).at(childconeiter->first);
                                const double locradius = levels_[level]->Cheb2Radius(childconeiter->first, locx);
                                std::complex<double> locfac;
                                _factorization(locradius, wavenumber_, locfac);

                                std::complex<double> new_locfac(1.0, 0.0); 

                                const double tmp_locfac_real = (locfac.real() * new_locfac.real() + locfac.imag() * new_locfac.imag()) / (new_locfac.real()*new_locfac.real() + new_locfac.imag()*new_locfac.imag());
                                const double tmp_locfac_imag = (-locfac.real() * new_locfac.imag() + locfac.imag() * new_locfac.real()) / (new_locfac.real()*new_locfac.real() + new_locfac.imag()*new_locfac.imag());

                                locfac = {tmp_locfac_real, tmp_locfac_imag};

                                std::vector<std::complex<double>> vals_begin;

                                if ((locrelconesegment >= levels_[level]->split_points_relconesegments_[comm_rank_]) && (locrelconesegment < levels_[level]->split_points_relconesegments_[comm_rank_+1])) {

                                    const long long coeffs_begin_id = (locrelconesegment - levels_[level]->split_points_relconesegments_[comm_rank_]) * P_; 
                                    vals_begin.insert(vals_begin.begin(), conesegments_prev.begin() + coeffs_begin_id, conesegments_prev.begin() + coeffs_begin_id + P_);

                                } else {

                                    vals_begin = coeffs.at(locrelconesegment);

                                }                             

                                const std::complex<double> tmp = IPSCHEME_.Interpolate(locx, locy, locz, &vals_begin[0]);

                                const double value_real = tmp.real() * locfac.real() - tmp.imag() * locfac.imag();
                                const double value_imag = tmp.imag() * locfac.real() + tmp.real() * locfac.imag();
                                
                                value[interpid] += std::complex<double>(value_real, value_imag);

                            }

                        }

                    }

                }

                for (int interpid = 0; interpid < IPSCHEME_.GetNInterpPoints(); interpid++) {

                    std::complex<double> fac;
                    _factorization(radii[interpid], wavenumber_, fac);
                    const double dabsfac = 1.0/(fac.real() * fac.real() + fac.imag() * fac.imag());

                    std::complex<double> new_fac(1.0, 0.0);    

                    const double tmp_fac_real = fac.real() * new_fac.real() + fac.imag() * new_fac.imag();
                    const double tmp_fac_imag = fac.real() * new_fac.imag() - fac.imag() * new_fac.real();

                    fac = std::complex<double>(tmp_fac_real, tmp_fac_imag);

                    const double value_real = (value[interpid].real() * fac.real() - value[interpid].imag() * fac.imag()) * dabsfac;
                    const double value_imag = (value[interpid].imag() * fac.real() + value[interpid].real() * fac.imag()) * dabsfac;
                    
                    conesegments_current[new_coeffs_begin_id+interpid] = {value_real, value_imag};

                }

                IPSCHEME_.GenerateInterpolant(&conesegments_current[new_coeffs_begin_id]);

            }

            }

        }

        void CommunicatePropagation(const int level, std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
                                    std::unordered_map<long long, long long> & num_idxs, std::unordered_map<long long, std::vector<long long>> & indexes,
                                    const std::vector<std::complex<double>>& conesegments_current) const 
        {

            coeffs.clear();
            num_idxs.clear();
            indexes.clear();

            for (int rank = 0; rank < comm_size_; rank++) {

                std::vector<std::complex<double>> coeffs_cone_segments_all;
                std::vector<std::complex<double>> coeffs_cone_segments_loc;

                std::vector<long long> num_idxs_loc;
                std::vector<long long> num_idxs_all;

                std::vector<long long> indexes_loc;
                std::vector<long long> indexes_all;

                std::vector<int> recv_counts_coeffs;
                std::vector<int> displs_coeffs;

                std::vector<int> recv_counts_num_idxs;
                std::vector<int> displs_num_idxs;

                if (rank == comm_rank_) {

                    for (const auto & cs : levels_[level]->cone_segments_propagation_not_in_rank_all_) {

                        const long long idx = (cs - levels_[level]->split_points_relconesegments_[rank]) * P_;

                        coeffs_cone_segments_all.insert(coeffs_cone_segments_all.end(), conesegments_current.begin() + idx, conesegments_current.begin() + idx + P_);
                        num_idxs_all.push_back(levels_[level]->num_elems_indexes_relconesegment_[idx / P_]);
                        indexes_all.insert(indexes_all.end(), levels_[level]->indexes_relconesegment_.begin() + idx, levels_[level]->indexes_relconesegment_.begin() + idx + P_);

                    }

                    recv_counts_coeffs = levels_[level]->recv_counts_propagation_;
                    displs_coeffs = levels_[level]->displs_propagation_;

                    recv_counts_num_idxs = levels_[level]->recv_counts_propagation_;
                    displs_num_idxs = levels_[level]->displs_propagation_;

                    for (int i = 0; i < comm_size_; i++) {
                        recv_counts_num_idxs[i] /= P_;
                        displs_num_idxs[i] /= P_;
                    }

                }

                int size_loc = levels_[level]->propagation_not_in_rank_[rank].size() * P_;
                
                coeffs_cone_segments_loc.resize(size_loc);
                num_idxs_loc.resize(size_loc / P_);
                indexes_loc.resize(size_loc);
                
                MPI_Scatterv(&coeffs_cone_segments_all[0], &recv_counts_coeffs[0], &displs_coeffs[0], MPI_DOUBLE_COMPLEX, &coeffs_cone_segments_loc[0], size_loc, MPI_DOUBLE_COMPLEX, rank, mpi_comm_);
                MPI_Scatterv(&indexes_all[0], &recv_counts_coeffs[0], &displs_coeffs[0], MPI_LONG_LONG, &indexes_loc[0], size_loc, MPI_LONG_LONG, rank, mpi_comm_);
                MPI_Scatterv(&num_idxs_all[0], &recv_counts_num_idxs[0], &displs_num_idxs[0], MPI_LONG_LONG, &num_idxs_loc[0], size_loc/P_, MPI_LONG_LONG, rank, mpi_comm_);

                for (long long i = 0; i < levels_[level]->propagation_not_in_rank_[rank].size(); i++) {

                    const long long cone_segment = levels_[level]->propagation_not_in_rank_[rank][i];

                    const std::vector<std::complex<double>> coefficients(coeffs_cone_segments_loc.begin() + i * P_, coeffs_cone_segments_loc.begin() + (i+1) * P_);
                    const std::vector<long long> indexes_vec(indexes_loc.begin() + i * P_, indexes_loc.begin() + (i+1) * P_);
                    const long long num_indx = num_idxs_loc[i];

                    coeffs[cone_segment] = coefficients;
                    num_idxs[cone_segment] = num_indx;
                    indexes[cone_segment] = indexes_vec;

                }

            }

        }

        // template<void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        // const double, const std::complex<double>, const std::complex<double>, std::complex<double>&)>
        // void Propagation(const int level, const std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
        //                  const std::unordered_map<long long, long long>& num_indxs, const std::unordered_map<long long, std::vector<long long>>& indexes,
        //                  const std::vector<std::complex<double>>& density,
        //                  std::vector<std::complex<double>>& conesegments_current, const std::vector<std::complex<double>>& conesegments_prev)
        void Propagation(
        int level,
        const std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
        const std::unordered_map<long long, long long>& num_indxs,
        const std::unordered_map<long long, std::vector<long long>>& indexes,
        const std::vector<std::complex<double>>& density,
        std::vector<std::complex<double>>& conesegments_current,
        const std::vector<std::complex<double>>& conesegments_prev,
        const std::function<void(
            double, double, double, double, double, double, double, double, double, std::complex<double>,
            std::complex<double>, std::complex<double>, std::complex<double>&)>& _kernel)
        {

            #pragma omp parallel
            {

            long long old_cocentered_morton_box = -1;
            std::vector<long long> morton_children;
            std::array<double, 3> center;
            
            #pragma omp for
            for (long long i = 0; i < levels_[level-1]->split_points_relconesegments_[comm_rank_+1]-levels_[level-1]->split_points_relconesegments_[comm_rank_]; i++) {

                const long long cocentered_morton_box = levels_[level-1]->relconesegment2cocenteredmortonboxid_[i];
                const long long nonrel_conesegment = levels_[level-1]->relconesegment2nonrelconesegment_[i];
                const double radius1 = levels_[level-1]->relconesegment2minradius_[i];
                const double radius2 = levels_[level-1]->relconesegment2maxradius_[i];

                Eigen::MatrixXd new_points = levels_[level-1]->interpolation_points_[nonrel_conesegment];

                const long long new_coeffs_begin_id = i * P_; 
                
                if (cocentered_morton_box != old_cocentered_morton_box) {

                    morton_children = GetMortonRelChildren(level-1, cocentered_morton_box);

                    old_cocentered_morton_box = cocentered_morton_box;

                    levels_[level-1]->GetBoxCenter(cocentered_morton_box, center[0], center[1], center[2]);

                } 

                Eigen::VectorXcd b(2*PS_II*PT_II);
                b.setZero();                

                for (int childiter = 0; childiter < morton_children.size(); childiter++) {

                    const long long morton_child = morton_children[childiter];
                    
                    std::array<long long, 2> data_box = levels_[level]->mortonbox2discretizationpoints_2_[morton_child];
                    
                    const long long points_begin = data_box[0];
                    const long long npoints = data_box[1];

                    for (int interpiter = 0; interpiter < PT_II*PS_II; interpiter++) {

                        const double x = center[0] + radius1 * new_points(0, interpiter);
                        const double y = center[1] + radius1 * new_points(1, interpiter);
                        const double z = center[2] + radius1 * new_points(2, interpiter);

                        if (USE_ADAPTIVITY && npoints <= MAX_ELEMS_LEAF) {

                            for (long long sourceiter = 0; sourceiter < npoints; sourceiter++) {

                                const double locx = x_[points_begin + sourceiter];
                                const double locy = y_[points_begin + sourceiter];
                                const double locz = z_[points_begin + sourceiter];

                                const double loc_normal_x = normal_x_[points_begin + sourceiter];
                                const double loc_normal_y = normal_y_[points_begin + sourceiter];
                                const double loc_normal_z = normal_z_[points_begin + sourceiter];

                                const std::complex<double> locdensity = density[points_begin + sourceiter];

                                std::complex<double> tmp_result;
                                _kernel(locx, locy, locz, x, y, z, loc_normal_x, loc_normal_y, loc_normal_z, coupling_parameter_, wavenumber_, locdensity, tmp_result);
                                    
                                b[interpiter] += tmp_result;

                            } 

                        } else {
                                
                            const long long locnonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_child, x, y, z);
                            const long long locrel_conesegment = levels_[level]->mortonboxnonrelcone2relcone_2_.at(morton_child).at(locnonrel_conesegment);

                            std::vector<std::complex<double>> vals_begin;
                            long long idx_num;
                            std::vector<long long> index;

                            if ((levels_[level]->split_points_relconesegments_[comm_rank_] <= locrel_conesegment) && (levels_[level]->split_points_relconesegments_[comm_rank_+1] > locrel_conesegment)) {

                                const long long coeffs_begin_id = (locrel_conesegment - levels_[level]->split_points_relconesegments_[comm_rank_]) * P_; 
                                vals_begin.insert(vals_begin.begin(), conesegments_prev.begin() + coeffs_begin_id, conesegments_prev.begin() + coeffs_begin_id + P_);
                                idx_num = levels_[level]->num_elems_indexes_relconesegment_[coeffs_begin_id / P_];
                                index.insert(index.begin(), levels_[level]->indexes_relconesegment_.begin() + coeffs_begin_id, levels_[level]->indexes_relconesegment_.begin() + coeffs_begin_id + P_);

                            } else {
                                
                                vals_begin = coeffs.at(locrel_conesegment);
                                idx_num = num_indxs.at(locrel_conesegment);
                                index = indexes.at(locrel_conesegment);
                                
                            }

                            for (long long ii = 0; ii < idx_num; ii++) {

                                const double locx = x_[index[ii]];
                                const double locy = y_[index[ii]];
                                const double locz = z_[index[ii]];

                                const double loc_normal_x = normal_x_[index[ii]];
                                const double loc_normal_y = normal_y_[index[ii]];
                                const double loc_normal_z = normal_z_[index[ii]];

                                const std::complex<double> locdensity = vals_begin[ii];

                                std::complex<double> tmp_result;
                                _kernel(locx, locy, locz, x, y, z, loc_normal_x, loc_normal_y, loc_normal_z, coupling_parameter_, wavenumber_, locdensity, tmp_result);
                        
                                b[interpiter] += tmp_result;

                            }

                        }

                    }

                }

                for (int childiter = 0; childiter < morton_children.size(); childiter++) {

                    const long long morton_child = morton_children[childiter];
                    
                    std::array<long long, 2> data_box = levels_[level]->mortonbox2discretizationpoints_2_[morton_child];
                    
                    const long long points_begin = data_box[0];
                    const long long npoints = data_box[1];

                    for (int interpiter = 0; interpiter < PT_II*PS_II; interpiter++) {

                        const double x = center[0] + radius2 * new_points(0, interpiter);
                        const double y = center[1] + radius2 * new_points(1, interpiter);
                        const double z = center[2] + radius2 * new_points(2, interpiter);

                        if (USE_ADAPTIVITY && npoints <= MAX_ELEMS_LEAF) {

                            for (long long sourceiter = 0; sourceiter < npoints; sourceiter++) {

                                const double locx = x_[points_begin + sourceiter];
                                const double locy = y_[points_begin + sourceiter];
                                const double locz = z_[points_begin + sourceiter];

                                const double loc_normal_x = normal_x_[points_begin + sourceiter];
                                const double loc_normal_y = normal_y_[points_begin + sourceiter];
                                const double loc_normal_z = normal_z_[points_begin + sourceiter];

                                const std::complex<double> locdensity = density[points_begin + sourceiter];

                                std::complex<double> tmp_result;
                                _kernel(locx, locy, locz, x, y, z, loc_normal_x, loc_normal_y, loc_normal_z, coupling_parameter_, wavenumber_, locdensity, tmp_result);
                                    
                                b[PS_II*PT_II + interpiter] += tmp_result;

                            } 

                        } else {
                                
                            const long long locnonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_child, x, y, z);
                            const long long locrel_conesegment = levels_[level]->mortonboxnonrelcone2relcone_2_.at(morton_child).at(locnonrel_conesegment);

                            std::vector<std::complex<double>> vals_begin;
                            long long idx_num;
                            std::vector<long long> index;

                            if ((levels_[level]->split_points_relconesegments_[comm_rank_] <= locrel_conesegment) && (levels_[level]->split_points_relconesegments_[comm_rank_+1] > locrel_conesegment)) {

                                const long long coeffs_begin_id = (locrel_conesegment - levels_[level]->split_points_relconesegments_[comm_rank_]) * P_; 
                                vals_begin.insert(vals_begin.begin(), conesegments_prev.begin() + coeffs_begin_id, conesegments_prev.begin() + coeffs_begin_id + P_);
                                idx_num = levels_[level]->num_elems_indexes_relconesegment_[coeffs_begin_id / P_];
                                index.insert(index.begin(), levels_[level]->indexes_relconesegment_.begin() + coeffs_begin_id, levels_[level]->indexes_relconesegment_.begin() + coeffs_begin_id + P_);

                            } else {
                                
                                vals_begin = coeffs.at(locrel_conesegment);
                                idx_num = num_indxs.at(locrel_conesegment);
                                index = indexes.at(locrel_conesegment);
                                
                            }

                            for (long long ii = 0; ii < idx_num; ii++) {

                                const double locx = x_[index[ii]];
                                const double locy = y_[index[ii]];
                                const double locz = z_[index[ii]];

                                const double loc_normal_x = normal_x_[index[ii]];
                                const double loc_normal_y = normal_y_[index[ii]];
                                const double loc_normal_z = normal_z_[index[ii]];

                                const std::complex<double> locdensity = vals_begin[ii];

                                std::complex<double> tmp_result;
                                _kernel(locx, locy, locz, x, y, z, loc_normal_x, loc_normal_y, loc_normal_z, coupling_parameter_, wavenumber_, locdensity, tmp_result);
                        
                                b[PS_II*PT_II + interpiter] += tmp_result;

                            }

                        }

                    }

                }

                const int num_idxs = levels_[level-1]->num_elems_indexes_relconesegment_[i];

                Eigen::MatrixXcd A(2*PS_II*PT_II, num_idxs);

                for (long long i1 = 0; i1 < PS_II; i1++) {
                    for (long long i2 = 0; i2 < PT_II; i2++) {

                        const double x = center[0] + radius1 * new_points(0, i1*PT_II+i2);
                        const double y = center[1] + radius1 * new_points(1, i1*PT_II+i2);
                        const double z = center[2] + radius1 * new_points(2, i1*PT_II+i2);

                        for (long long jj = 0; jj < num_idxs; jj++) {

                            const long long idx = levels_[level-1]->indexes_relconesegment_[new_coeffs_begin_id + jj];

                            std::complex<double> kernel;   

                            const double locx = x_[idx];
                            const double locy = y_[idx];
                            const double locz = z_[idx];

                            const double loc_normal_x = normal_x_[idx];
                            const double loc_normal_y = normal_y_[idx];
                            const double loc_normal_z = normal_z_[idx];

                            const std::complex<double> locdensity = {1.0, 0.0};

                            _kernel(locx, locy, locz,
                                    x, y, z,
                                    loc_normal_x, loc_normal_y, loc_normal_z,
                                    coupling_parameter_, wavenumber_,
                                    locdensity,
                                    kernel);

                            A(i1 * PT_II + i2, jj) = kernel;

                        }

                    }
                } 

                for (long long i1 = 0; i1 < PS_II; i1++) {
                    for (long long i2 = 0; i2 < PT_II; i2++) {

                        const double x = center[0] + radius2 * new_points(0, i1*PT_II+i2);
                        const double y = center[1] + radius2 * new_points(1, i1*PT_II+i2);
                        const double z = center[2] + radius2 * new_points(2, i1*PT_II+i2);

                        for (long long jj = 0; jj < num_idxs; jj++) {

                            const long long idx = levels_[level-1]->indexes_relconesegment_[new_coeffs_begin_id + jj];

                            std::complex<double> kernel;   

                            const double locx = x_[idx];
                            const double locy = y_[idx];
                            const double locz = z_[idx];

                            const double loc_normal_x = normal_x_[idx];
                            const double loc_normal_y = normal_y_[idx];
                            const double loc_normal_z = normal_z_[idx];

                            const std::complex<double> locdensity = {1.0, 0.0};

                            _kernel(locx, locy, locz,
                                    x, y, z,
                                    loc_normal_x, loc_normal_y, loc_normal_z,
                                    coupling_parameter_, wavenumber_,
                                    locdensity,
                                    kernel);

                            A(PS_II*PT_II + i1 * PT_II + i2, jj) = kernel;

                        }

                    }
                } 

                Eigen::VectorXcd coeffs = A.completeOrthogonalDecomposition().solve(b);
                
                for (long long ii = 0; ii < num_idxs; ii++) {

                    conesegments_current[new_coeffs_begin_id + ii] = coeffs[ii];

                }
                
            }

            }

        }      

    public:

        BoxTree()
        {

        }

        ~BoxTree() 
        {

            for (int iter = 0; iter < nlevels_; iter++) {

                delete levels_[iter];

            }

        }

        void InitializeObject(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, 
                              const std::vector<double>& normal_x, const std::vector<double>& normal_y, const std::vector<double>& normal_z,
                              const std::vector<long long>& split_points, const std::vector<int>& recv_counts, const std::vector<int>& displs,
                              const std::complex<double> coupling_parameter, std::vector<long long>& sorting, int nlevels, std::complex<double> wavenumber,
                              const MPI_Comm& mpi_comm)
                              
        {
            mpi_comm_ = mpi_comm;

            x_ = x;
            y_ = y;
            z_ = z;
            
            normal_x_ = normal_x;
            normal_y_ = normal_y;
            normal_z_ = normal_z;
            
            nlevels_ = nlevels;
            wavenumber_ = wavenumber;

            split_points_orig_ = split_points;
            recv_counts_orig_ = recv_counts;
            displs_orig_ = displs;

            N_loc_orig_ = x.size();

            coupling_parameter_ = coupling_parameter;

            if (USE_HIGH_ORDER) {
                if (TOTAL_COLS != -1) {
                    P_ = std::min(TOTAL_COLS, 2*PS_II*PT_II);
                } else {
                    P_ = 2*PS_II*PT_II;
                }
            } else {
                P_ = PS * PT * PT;
            }

            Initialize();

            sorting = sorting_;

            x = x_;
            y = y_;
            z = z_;

        }
        
        long long get_box(const long long i)
        {

            return levels_.back()->Point2Morton(x_[i], y_[i], z_[i]);

        }   

        void get_mortonidofrelboxes(std::unordered_set<long long> & ids) 
        {

            ids = levels_.back()->mortonidofrelboxes_;

        }

        std::vector<long long> get_neighbours_box(const long long i) 
        {

            return levels_.back()->GetNeighboursAll(i);

        }

        void set_precomputations_data(const std::unordered_map<long long, std::unordered_set<long long>>& data)
        {

            precomputations_data_ = data;

        }

        void set_n_pts_per_patch(const int n1, const int n2) {

            n_pts_per_patch_.resize(2);

            n_pts_per_patch_[0] = n1;
            n_pts_per_patch_[1] = n2;

        }

        // template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        // const double, const std::complex<double>, const std::complex<double>, std::complex<double>&)>
        // void InitializeIndexes() {

        void InitializeIndexes(
        const std::function<void(
        double, double, double, double, double, double, double, double, double, std::complex<double>,
        std::complex<double>, std::complex<double>, std::complex<double>&)>& _kernel) {

            const int D = nlevels_ - 1; 

            std::unordered_map<long long, std::vector<long long>> tmp_indx;                

            // LevelDIndexes<_kernel>(); 
            LevelDIndexes(_kernel);

            CommunicateIndexes(D, tmp_indx);

            for (int level = D; level >= 2; level--) {

                if (level > 2) {

                    LeveldIndexes(level, tmp_indx, _kernel);

                    // LeveldIndexes<_kernel>(level, tmp_indx);

                    if (level > 3) {

                        CommunicateIndexes(level-1, tmp_indx);

                    }

                }

            }

        }


        void Solve(
            std::vector<std::complex<double>>& density,
            const std::function<void(
                double, double, double, double, double, double, double, double, double, std::complex<double>,
                std::complex<double>, std::complex<double>, std::complex<double>&)>& kernel,
            const std::function<void(
                double, std::complex<double>, std::complex<double>&)>& factorization)
        {
            if (USE_HIGH_ORDER) {
                SolveHO(density, kernel);
            } else {
                SolveOrig(density, kernel, factorization);
            }
        }

        // template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        // const double, const std::complex<double>, const std::complex<double>, std::complex<double>&),
        // void _factorization(const double, const std::complex<double>, std::complex<double>&)>
        // void Solve(std::vector<std::complex<double>>& density) {

        //     if (USE_HIGH_ORDER) {

        //         SolveHO<_kernel>(density);

        //     } else {

        //         SolveOrig<_kernel, _factorization>(density);

        //     }

        // }

                // template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        // const double, const std::complex<double>, const std::complex<double>, std::complex<double>&),
        // void _factorization(const double, const std::complex<double>, std::complex<double>&)>
        // void SolveOrig(std::vector<std::complex<double>>& density) {

        void SolveOrig(
    std::vector<std::complex<double>>& density,
    const std::function<void(
        double, double, double, double, double, double, double, double, double, std::complex<double>,
        std::complex<double>, std::complex<double>, std::complex<double>&)>& _kernel,
        const std::function<void(
        double, std::complex<double>, std::complex<double>&)>& _factorization) { 



            if (nlevels_ < 3) 
                throw std::invalid_argument("IFGF accelerator requires at least 3 levels.");
            
            long long maxncoeffs = 0;

            for (int level = 2; level < nlevels_; level++) {

                maxncoeffs = std::max<long long>(maxncoeffs, levels_[level]->relconesegment2nonrelconesegment_.size());

            }

            std::vector<std::complex<double>> conesegments_current = std::vector<std::complex<double>>(maxncoeffs * P_, {0.0, 0.0});
            std::vector<std::complex<double>> conesegments_prev = std::vector<std::complex<double>>(maxncoeffs * P_, {0.0, 0.0});          
            
            const int D = nlevels_ - 1;   

            std::unordered_map<long long, std::vector<std::complex<double>>> tmpinterpolationcoeffs;         
            std::unordered_map<long long, std::vector<std::complex<double>>> tmppropagationcoeffs;         

            ZeroSolution();

            // SingularInteractions<_kernel>(density);
            SingularInteractions(density, _kernel);

            // LevelDEvaluations<_kernel, _factorization>(density, conesegments_current, conesegments_prev);
            LevelDEvaluations(density, conesegments_current, conesegments_prev, _kernel, _factorization);
            
            CommunicatePropagation(D, tmppropagationcoeffs, conesegments_current);            

            for (int level = D; level >= 2; level--) {

                CommunicateInterpolation(level, tmpinterpolationcoeffs, conesegments_prev);

                if (level > 2) {
                    
                    // Propagation<_kernel, _factorization>(level, tmppropagationcoeffs, density, conesegments_current, conesegments_prev);
                         
                    Propagation(level, tmppropagationcoeffs, density, conesegments_current, conesegments_prev, _kernel, _factorization);

                    if (level > 3) {

                        CommunicatePropagation(level-1, tmppropagationcoeffs, conesegments_current);

                    }

                }             

                // Interpolation<_kernel, _factorization>(level, tmpinterpolationcoeffs, density, conesegments_prev);
                Interpolation(level, tmpinterpolationcoeffs, density, conesegments_prev, _kernel, _factorization);

                SwapCoefficients(conesegments_current, conesegments_prev);

            }


            std::vector<std::complex<double>>().swap(conesegments_current);
            std::vector<std::complex<double>>().swap(conesegments_prev);

            tmppropagationcoeffs.clear();
            tmpinterpolationcoeffs.clear();

            const long long N = x_.size();

            std::vector<std::complex<double>> solution_all(N);

            MPI_Allgatherv(&solution_[0], N_loc_orig_, MPI_DOUBLE_COMPLEX, &solution_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE_COMPLEX, mpi_comm_);
            
            density = solution_all;            

        }

                // template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        // const double, const std::complex<double>, const std::complex<double>, std::complex<double>&)>
        // void SolveHO(std::vector<std::complex<double>>& density) {

        void SolveHO(
    std::vector<std::complex<double>>& density,
    const std::function<void(
        double, double, double, double, double, double, double, double, double, std::complex<double>,
        std::complex<double>, std::complex<double>, std::complex<double>&)>& _kernel)
        {

            if (nlevels_ < 3) 
                throw std::invalid_argument("IFGF accelerator requires at least 3 levels.");
            
            long long maxncoeffs = 0;

            for (int level = 2; level < nlevels_; level++) {

                maxncoeffs = std::max<long long>(maxncoeffs, levels_[level]->relconesegment2nonrelconesegment_.size());

            }

            std::vector<std::complex<double>> conesegments_current = std::vector<std::complex<double>>(maxncoeffs * P_, {0.0, 0.0});
            std::vector<std::complex<double>> conesegments_prev = std::vector<std::complex<double>>(maxncoeffs * P_, {0.0, 0.0});

            const int D = nlevels_ - 1;   

            std::unordered_map<long long, std::vector<std::complex<double>>> tmpinterpolationcoeffs;         
            std::unordered_map<long long, std::vector<std::complex<double>>> tmppropagationcoeffs;   

            std::unordered_map<long long, long long> tmpinterpolationnum_indx;
            std::unordered_map<long long, long long> tmppropagationnum_indx;

            std::unordered_map<long long, std::vector<long long>> tmpinterpolationindexes;
            std::unordered_map<long long, std::vector<long long>> tmppropagationindexes;      
            
            ZeroSolution();

            // SingularInteractions<_kernel>(density);
               SingularInteractions(density, _kernel);

            // LevelDEvaluations<_kernel>(density, conesegments_current, conesegments_prev);
            LevelDEvaluations(density, conesegments_current, conesegments_prev, _kernel);

            CommunicatePropagation(D, tmppropagationcoeffs, tmppropagationnum_indx, tmppropagationindexes, conesegments_current);    

            for (int level = D; level >= 2; level--) {

                CommunicateInterpolation(level, tmpinterpolationcoeffs, tmpinterpolationnum_indx, tmpinterpolationindexes, conesegments_prev);

                if (level > 2) {
                    
                    // Propagation<_kernel>(level, tmppropagationcoeffs, tmppropagationnum_indx, tmppropagationindexes, density, conesegments_current, conesegments_prev);
                    Propagation(level, tmppropagationcoeffs, tmppropagationnum_indx, tmppropagationindexes, density, conesegments_current, conesegments_prev, _kernel);

                    if (level > 3) {

                        CommunicatePropagation(level-1, tmppropagationcoeffs, tmppropagationnum_indx, tmppropagationindexes, conesegments_current);          

                    }

                }             

                // Interpolation<_kernel>(level, tmpinterpolationcoeffs, tmpinterpolationnum_indx, tmpinterpolationindexes, density, conesegments_prev);
                  Interpolation(level, tmpinterpolationcoeffs, tmpinterpolationnum_indx, tmpinterpolationindexes, density, conesegments_prev, _kernel);

                SwapCoefficients(conesegments_current, conesegments_prev);

            }  

            const long long N = x_.size();

            std::vector<std::complex<double>> solution_all(N);

            MPI_Allgatherv(&solution_[0], N_loc_orig_, MPI_DOUBLE_COMPLEX, &solution_all[0], &recv_counts_orig_[0], &displs_orig_[0], MPI_DOUBLE_COMPLEX, mpi_comm_);
                        
            density = solution_all;              

        }

};

#endif