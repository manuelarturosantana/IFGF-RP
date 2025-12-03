#ifndef BOXTREE_H
#define BOXTREE_H

#include <vector>
#include <limits>
#include <map>
#include <set>
#include <iostream>
#include <numeric>
#include <algorithm>

#include "mkl.h"
#include <omp.h>
#include "mpi.h"

#include "Level.h"

#include "Interpolator.h"

#include "solver2.h"

#include "global.h"

#include <fstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <cstring>

class Level;

template <int PS, int PT>
class BoxTree;

template <int PS, int PT>
class BoxTree 
{

    private:

        int comm_rank_;
        int comm_size_;
        MPI_Comm mpi_comm_;

        std::vector<int> N_PTS_PER_PATCH = std::vector<int>(2);

        std::vector<long long> split_points_orig_;
        std::vector<MPI_Count> recv_counts_orig_;
        std::vector<MPI_Aint> displs_orig_;

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

        std::vector<std::unordered_set<long long>> points_not_in_rank_;
        std::unordered_map<long long, std::array<double, 7>> data_not_in_rank_;

        // These parameters need to be passed in from Solver
        int nlevels_;
        std::complex<double> wavenumber_;
        std::complex<double> coupling_parameter_; 
        double BBSIZEOFFSET;
        bool USE_ADAPTIVITY;
        bool USE_ACCELERATOR;
        long long MAX_ELEMS_LEAF;


     

        std::unordered_map<long long, std::unordered_set<long long>> precomputations_data_;

        std::vector<Level*> levels_;

        std::vector<std::complex<double>> solution_;
        
        int P_;

        long long max_elems_level_D_;

        Interpolator<PS, PT> IPSCHEME_;

    private:

        void Initialize() 
        {

            InitializeMPI();

            MPI_Barrier(mpi_comm_);

            InitializeBoxesAndLevels();

            MPI_Barrier(mpi_comm_);

            InitializeRelevantConeSegments();

            MPI_Barrier(mpi_comm_);

            InitializePointsNotInRank();

            MPI_Barrier(mpi_comm_);

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

        void ComputeNLevels(const std::array<double, 3>& min, const double boxsize_x, const double boxsize_y, const double boxsize_z)
        {

            int nlevels = 0;

            bool split = true;

            while (split) {

                nlevels++;

                const long long nboxes_axis = 1LL << (nlevels - 1);
                
                std::unordered_map<long long, long long> morton_code_counts_loc;                

                for (long long i = 0; i < N_loc_orig_; i++) {

                    const long long morton_box = Level::Point2Morton(x_[i], y_[i], z_[i], min[0], min[1], min[2], boxsize_x / nboxes_axis, boxsize_y / nboxes_axis, boxsize_z / nboxes_axis, nlevels - 1);
 
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

                MPI_Allgather(&size_loc, 1, MPI_COUNT, &recv_counts[0], 1, MPI_COUNT, mpi_comm_);

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

            nlevels_ = nlevels-1;

            if (comm_rank_ == 0) {

                std::cout << "Number of levels = " << nlevels_ << "\n";

            }

        }

        void CountingSortDigit(std::vector<long long>& data, std::vector<long long>& order, long long exp) 
        {

            long long n = data.size();

            std::vector<long long> output(n);
            std::vector<long long> output_order(n);
            std::vector<long long> count(10, 0);
            
            for (long long i = 0; i < n; i++) {

                count[(data[i] / exp) % 10]++;

            }
        
            for (int i = 1; i < 10; i++) {

                count[i] += count[i - 1];

            }
        
            for (long long i = n-1; i >= 0; i--) {

                output[count[(data[i] / exp) % 10] - 1] = data[i];
                output_order[count[(data[i] / exp) % 10] - 1] = order[i];
                count[(data[i] / exp) % 10]--;

            }
        
            std::copy(output.begin(), output.end(), data.begin());
            std::copy(output_order.begin(), output_order.end(), order.begin());

        }
        
        void ParallelRadixSort(std::vector<long long>& local_data, std::vector<long long>& local_order)
        {
            
            long long local_max = *std::max_element(local_data.begin(), local_data.end());
            long long global_max;
            MPI_Allreduce(&local_max, &global_max, 1, MPI_LONG_LONG, MPI_MAX, mpi_comm_);
        
            for (long long exp = 1; global_max / exp > 0; exp *= 10) {

                CountingSortDigit(local_data, local_order, exp); 
                
                std::vector<long long> local_counts(10, 0);  

                for (long long val : local_data) {

                    local_counts[(val / exp) % 10]++;

                }    
                
                std::vector<long long> global_counts(10, 0);
                MPI_Allreduce(local_counts.data(), global_counts.data(), 10, MPI_LONG_LONG, MPI_SUM, mpi_comm_);

                std::vector<long long> partial_counts(10, 0);
                MPI_Exscan(local_counts.data(), partial_counts.data(), 10, MPI_LONG_LONG, MPI_SUM, mpi_comm_);

                std::vector<long long> global_displs(11, 0);
                for (int i = 1; i < 11; i++) {
                    global_displs[i] = global_displs[i-1] + global_counts[i-1];
                }

                std::vector<MPI_Count> send_counts(comm_size_, 0);
                std::vector<MPI_Aint> send_displs(comm_size_, 0);                

                for (int i = 0; i < comm_size_; ++i) {

                    long long pos_low = split_points_orig_[i];
                    long long pos_up = split_points_orig_[i+1]-1;
                    
                    auto it_low = std::upper_bound(global_displs.begin(), global_displs.end(), pos_low);
                    int digit_low = static_cast<int>(std::distance(global_displs.begin(), it_low)) - 1;
                    if (digit_low < 0) { 
                        digit_low = 0;            
                    }
                    if (digit_low >= 10) {                     
                        digit_low = 9; 
                    }

                    auto it_up = std::upper_bound(global_displs.begin(), global_displs.end(), pos_up);
                    int digit_up = static_cast<int>(std::distance(global_displs.begin(), it_up)) - 1;
                    if (digit_up < 0) { 
                        digit_up = 0;            
                    }
                    if (digit_up >= 10) {                     
                        digit_up = 9; 
                    }

                    for (int j = digit_low; j <= digit_up; j++) {

                        if (local_counts[j] == 0) {
                            continue;
                        }

                        long long needed_count;

                        if (digit_low == digit_up) {
                            needed_count = pos_up + 1 - pos_low;
                        } else if (j == digit_low) {
                            needed_count = global_displs[j+1] - pos_low;
                        } else if (j == digit_up) {
                            needed_count = pos_up + 1 - global_displs[j];
                        } else {
                            needed_count = global_displs[j+1] - global_displs[j];
                        }

                        if (partial_counts[j] >= needed_count) {
                            partial_counts[j] -= needed_count;
                        } else {
                            long long total_present = local_counts[j];
                            long long already_sent = partial_counts[j];
                            long long need_to_be_sent = needed_count - already_sent;
                            long long total_sent = std::min(total_present, need_to_be_sent); 
                            send_counts[i] += total_sent;
                            local_counts[j] -= total_sent; 
                            partial_counts[j] = 0;
                        }                           

                    }

                    if (i > 0) {

                        send_displs[i] = send_displs[i-1] + send_counts[i-1];

                    }

                }

                std::vector<MPI_Count> recv_counts(comm_size_, 0);
                std::vector<MPI_Aint> recv_displs(comm_size_, 0);

                MPI_Alltoall(send_counts.data(), 1, MPI_COUNT, recv_counts.data(), 1, MPI_COUNT, mpi_comm_);
                
                for (int i = 1; i < comm_size_; ++i) {

                    recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];

                }

                std::vector<long long> new_local_data(N_loc_orig_);
                std::vector<long long> new_local_order(N_loc_orig_);

                MPI_Alltoallv_c(local_data.data(), send_counts.data(), send_displs.data(), MPI_LONG_LONG,
                                new_local_data.data(), recv_counts.data(), recv_displs.data(), MPI_LONG_LONG,
                                mpi_comm_);
                MPI_Alltoallv_c(local_order.data(), send_counts.data(), send_displs.data(), MPI_LONG_LONG,
                                new_local_order.data(), recv_counts.data(), recv_displs.data(), MPI_LONG_LONG,
                                mpi_comm_);

                local_data.swap(new_local_data);
                local_order.swap(new_local_order);

                CountingSortDigit(local_data, local_order, exp);

            }           

        }

        void SortBox(const std::array<double, 3>& min, const double boxsize_x, const double boxsize_y, const double boxsize_z) 
        {

            std::vector<long long> morton_code(N_loc_orig_);

            #pragma omp parallel for
            for (long long i = 0; i < N_loc_orig_; i++) {

                const long long morton = Level::Point2Morton(x_[i], y_[i], z_[i], min[0], min[1], min[2], boxsize_x, boxsize_y, boxsize_z, nlevels_-1);
                morton_code[i] = morton;

            } 

            ParallelRadixSort(morton_code, sorting_);

            std::vector<long long>().swap(morton_code);

            std::vector<MPI_Count> send_counts(comm_size_, 0);

            for (long long i = 0; i < N_loc_orig_; i++) {

                auto it = std::upper_bound(split_points_orig_.begin(), split_points_orig_.end(), sorting_[i]);

                int rank = std::clamp<int>(static_cast<int>(std::distance(split_points_orig_.begin(), it)) - 1, 0, comm_size_-1);

                send_counts[rank]++;

            }

            std::vector<MPI_Count> recv_counts(comm_size_);

            MPI_Alltoall(send_counts.data(), 1, MPI_COUNT, recv_counts.data(), 1, MPI_COUNT, mpi_comm_);

            std::vector<MPI_Aint> send_displs(comm_size_), recv_displs(comm_size_);

            std::exclusive_scan(send_counts.begin(), send_counts.end(), send_displs.begin(), 0);
            std::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_displs.begin(), 0);

            std::vector<long long> send_idx(N_loc_orig_);
            std::vector<MPI_Aint> local_offsets = send_displs;

            std::vector<long long> indexes_map(N_loc_orig_);

            for (long long i = 0; i < N_loc_orig_; ++i) {

                auto it = std::upper_bound(split_points_orig_.begin(), split_points_orig_.end(), sorting_[i]);
                
                int rank = std::clamp<int>(static_cast<int>(std::distance(split_points_orig_.begin(), it)) - 1, 0, comm_size_ - 1);

                long long pos = local_offsets[rank];
                send_idx[pos] = sorting_[i];
                local_offsets[rank]++;

                indexes_map[i] = pos;

            }

            std::vector<double> x_send(N_loc_orig_), y_send(N_loc_orig_), z_send(N_loc_orig_);
            std::vector<double> nx_send(N_loc_orig_), ny_send(N_loc_orig_), nz_send(N_loc_orig_);
            std::vector<long long> sorting_send(N_loc_orig_);

            MPI_Alltoallv_c(send_idx.data(), send_counts.data(), send_displs.data(), MPI_LONG_LONG,
                            sorting_send.data(), recv_counts.data(), recv_displs.data(), MPI_LONG_LONG,
                            mpi_comm_);

            std::vector<long long>().swap(send_idx);

            #pragma omp parallel for
            for (long long i = 0; i < N_loc_orig_; ++i) {

                const long long j = sorting_send[i] - point_low_;

                x_send[i]  = x_[j];
                y_send[i]  = y_[j];
                z_send[i]  = z_[j];

                nx_send[i] = normal_x_[j];
                ny_send[i] = normal_y_[j];
                nz_send[i] = normal_z_[j];

            }

            std::vector<long long>().swap(sorting_send);

            MPI_Alltoallv_c(x_send.data(), recv_counts.data(), recv_displs.data(), MPI_DOUBLE,
                            x_.data(), send_counts.data(), send_displs.data(), MPI_DOUBLE,
                            mpi_comm_);
            MPI_Alltoallv_c(y_send.data(), recv_counts.data(), recv_displs.data(), MPI_DOUBLE,
                            y_.data(), send_counts.data(), send_displs.data(), MPI_DOUBLE,
                            mpi_comm_);
            MPI_Alltoallv_c(z_send.data(), recv_counts.data(), recv_displs.data(), MPI_DOUBLE,
                            z_.data(), send_counts.data(), send_displs.data(), MPI_DOUBLE,
                            mpi_comm_);
            MPI_Alltoallv_c(nx_send.data(), recv_counts.data(), recv_displs.data(), MPI_DOUBLE,
                            normal_x_.data(), send_counts.data(), send_displs.data(), MPI_DOUBLE,
                            mpi_comm_);
            MPI_Alltoallv_c(ny_send.data(), recv_counts.data(), recv_displs.data(), MPI_DOUBLE,
                            normal_y_.data(), send_counts.data(), send_displs.data(), MPI_DOUBLE,
                            mpi_comm_);
            MPI_Alltoallv_c(nz_send.data(), recv_counts.data(), recv_displs.data(), MPI_DOUBLE,
                            normal_z_.data(), send_counts.data(), send_displs.data(), MPI_DOUBLE,
                            mpi_comm_);

            #pragma omp parallel for
            for (long long i = 0; i < N_loc_orig_; i++) {

                long long pos = indexes_map[i];

                x_send[i] = x_[pos];
                y_send[i] = y_[pos];
                z_send[i] = z_[pos];

                nx_send[i] = normal_x_[pos];
                ny_send[i] = normal_y_[pos];
                nz_send[i] = normal_z_[pos];
                    
            }

            x_ = std::move(x_send);
            y_ = std::move(y_send);
            z_ = std::move(z_send);

            normal_x_ = std::move(nx_send);
            normal_y_ = std::move(ny_send);
            normal_z_ = std::move(nz_send);

            std::vector<double>().swap(x_send);
            std::vector<double>().swap(y_send);
            std::vector<double>().swap(z_send);

            std::vector<double>().swap(nx_send);
            std::vector<double>().swap(ny_send);
            std::vector<double>().swap(nz_send);
            
            std::vector<long long>().swap(indexes_map); 

        }

        void InitializeLevelDBoxesData() 
        {

            long long old_morton = -1;
            long long npoints;

            for (long long i = 0; i < N_loc_orig_; i++) {

                const long long idx = i + point_low_;

                const long long morton = levels_.back()->Point2Morton(x_[i], y_[i], z_[i]);

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

            std::vector<long long> relevantmorton_loc(levels_.back()->mortonidofrelboxes_.begin(), levels_.back()->mortonidofrelboxes_.end());
            MPI_Count numrelboxes_loc = relevantmorton_loc.size();

            std::vector<MPI_Count> numrelboxes_all(comm_size_);
            MPI_Allgather(&numrelboxes_loc, 1, MPI_COUNT, &numrelboxes_all[0], 1, MPI_COUNT, mpi_comm_);

            std::vector<MPI_Aint> displs(comm_size_);
            MPI_Count numrelboxes_total = 0;

            for (int i = 0; i < comm_size_; i++) {

                displs[i] = numrelboxes_total;
                numrelboxes_total += numrelboxes_all[i];

            }
            
            std::vector<long long> mortonidofrelboxes_all(static_cast<size_t>(numrelboxes_total));
            MPI_Allgatherv_c(&relevantmorton_loc[0], numrelboxes_loc, MPI_LONG_LONG, &mortonidofrelboxes_all[0], &numrelboxes_all[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
       
            std::vector<long long> relevantmortonown(levels_.back()->mortonidofrelboxes_.begin(), levels_.back()->mortonidofrelboxes_.end());
            std::sort(relevantmortonown.begin(), relevantmortonown.end());

            std::vector<long long> relevantmortonothers;
            relevantmortonothers.reserve(mortonidofrelboxes_all.size());

            for (long long morton_box : mortonidofrelboxes_all) {

                bool is_own = std::binary_search(relevantmortonown.begin(), relevantmortonown.end(), morton_box);

                if (!is_own)
                    relevantmortonothers.push_back(morton_box);

            }

            std::sort(relevantmortonothers.begin(), relevantmortonothers.end());
            relevantmortonothers.erase(std::unique(relevantmortonothers.begin(), relevantmortonothers.end()), relevantmortonothers.end());

            std::vector<long long>().swap(relevantmorton_loc);
            std::vector<long long>().swap(mortonidofrelboxes_all);

            std::vector<long long> relevantparentmortonown;
            std::vector<long long> relevantparentmortonothers;

            relevantparentmortonown.reserve(relevantmortonown.size());
            relevantparentmortonothers.reserve(relevantmortonothers.size());

            for (int i = nlevels_-1; i >= 2; i--) {

                if (i != nlevels_-1) {

                    levels_[i]->mortonidofrelboxes_.insert(relevantmortonown.begin(), relevantmortonown.end());

                }

                for (const auto & morton_box : relevantmortonown) {

                    std::vector<long long> neighbors_and_cousins = levels_[i]->GetAllNeighboursAndCousins(morton_box);

                    for (long long nc_morton_box : neighbors_and_cousins) {

                        bool in_own = std::binary_search(relevantmortonown.begin(), relevantmortonown.end(), nc_morton_box);
                        bool in_other = std::binary_search(relevantmortonothers.begin(), relevantmortonothers.end(), nc_morton_box);

                        if (in_own || in_other) {
                            levels_[i]->mortonidofrelboxes_2_.insert(nc_morton_box);
                        }

                    }

                    relevantparentmortonown.push_back(static_cast<long long>(morton_box/8));

                }

                for (const auto & morton_box : relevantmortonothers) {

                    const long long parent_morton = static_cast<long long>(morton_box/8);

                    if (!std::binary_search(relevantparentmortonown.begin(), relevantparentmortonown.end(), parent_morton)) {

                        relevantparentmortonothers.push_back(parent_morton);

                    }

                }

                std::sort(relevantparentmortonown.begin(), relevantparentmortonown.end());
                relevantparentmortonown.erase(std::unique(relevantparentmortonown.begin(), relevantparentmortonown.end()), relevantparentmortonown.end());

                std::sort(relevantparentmortonothers.begin(), relevantparentmortonothers.end());
                relevantparentmortonothers.erase(std::unique(relevantparentmortonothers.begin(), relevantparentmortonothers.end()), relevantparentmortonothers.end());

                std::swap(relevantmortonown, relevantparentmortonown);
                std::swap(relevantmortonothers, relevantparentmortonothers);

                relevantparentmortonown.clear();
                relevantparentmortonothers.clear();

            }   
            
        }

        void InitializeBoxesDataAllLevels() 
        {

            std::vector<long long> morton_box_loc;
            std::vector<long long> position_loc;
            std::vector<long long> size_loc;

            morton_box_loc.reserve(levels_.back()->mortonbox2discretizationpoints_.size());
            position_loc.reserve(levels_.back()->mortonbox2discretizationpoints_.size());
            size_loc.reserve(levels_.back()->mortonbox2discretizationpoints_.size());

            for (const auto & [morton_box, values] : levels_.back()->mortonbox2discretizationpoints_) {

                morton_box_loc.push_back(morton_box);
                position_loc.push_back(values[0]);
                size_loc.push_back(values[1]);

            }

            MPI_Count total_size_loc = morton_box_loc.size();

            std::vector<MPI_Count> recv_counts(comm_size_);
            std::vector<MPI_Aint> displs(comm_size_);

            MPI_Allgather(&total_size_loc, 1, MPI_COUNT, &recv_counts[0], 1, MPI_COUNT, mpi_comm_);

            MPI_Count total_size = 0;

            for (int i = 0; i < comm_size_; i++) {

                displs[i] = total_size;
                total_size += recv_counts[i];

            }

            std::vector<long long> morton_box_all(static_cast<size_t>(total_size));
            std::vector<long long> position_all(static_cast<size_t>(total_size));
            std::vector<long long> size_all(static_cast<size_t>(total_size));

            MPI_Allgatherv_c(&morton_box_loc[0], total_size_loc, MPI_LONG_LONG, &morton_box_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
            MPI_Allgatherv_c(&position_loc[0], total_size_loc, MPI_LONG_LONG, &position_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
            MPI_Allgatherv_c(&size_loc[0], total_size_loc, MPI_LONG_LONG, &size_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);

            std::vector<long long>().swap(morton_box_loc);
            std::vector<long long>().swap(position_loc);
            std::vector<long long>().swap(size_loc);
            
            std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_all;

            max_elems_level_D_ = 0;

            for (long long i = 0; i < total_size; i++) {

                const long long morton_box = morton_box_all[i];
                const long long position = position_all[i];
                const long long size = size_all[i];

                auto [it, inserted] = mortonbox2discretizationpoints_all.try_emplace(morton_box, std::array<long long, 2>{position, size});

                if (!inserted) {

                    auto & entry = it->second;

                    entry[0] = std::min(entry[0], position);
                    entry[1] += size;

                }

                max_elems_level_D_ = std::max(max_elems_level_D_, size);

            }             

            std::vector<long long>().swap(morton_box_all);
            std::vector<long long>().swap(position_all);
            std::vector<long long>().swap(size_all); 

            std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_loc = levels_.back()->mortonbox2discretizationpoints_;
            
            std::unordered_map<long long, std::array<long long, 2>> new_mortonbox2discretizationpoints_all;
            std::unordered_map<long long, std::array<long long, 2>> new_mortonbox2discretizationpoints_loc;

            for (int i = nlevels_-1; i >= 2; i--) {

                if (i != nlevels_-1) {

                    levels_[i]->mortonbox2discretizationpoints_ = mortonbox2discretizationpoints_loc;

                }

                for (const auto & morton_box : levels_[i]->mortonidofrelboxes_2_) {

                    levels_[i]->mortonbox2discretizationpoints_2_[morton_box] = mortonbox2discretizationpoints_all[morton_box];

                }

                new_mortonbox2discretizationpoints_all.clear();
                new_mortonbox2discretizationpoints_all.reserve(mortonbox2discretizationpoints_all.size() / 8);
                
                for (const auto [morton_box, values] : mortonbox2discretizationpoints_all) {

                    const long long new_morton_box = static_cast<long long>(morton_box/8);

                    auto [it, inserted] = new_mortonbox2discretizationpoints_all.try_emplace(new_morton_box, values);

                    if (!inserted) {

                        auto& entry = it->second;

                        entry[0] = std::min(entry[0], values[0]);
                        entry[1] += values[1];

                    }

                }

                mortonbox2discretizationpoints_all = new_mortonbox2discretizationpoints_all;

                new_mortonbox2discretizationpoints_loc.clear();
                new_mortonbox2discretizationpoints_loc.reserve(mortonbox2discretizationpoints_loc.size() / 8);

                for (const auto [morton_box, values] : mortonbox2discretizationpoints_loc) {

                    const long long new_morton_box = static_cast<long long>(morton_box/8);

                    auto [it, inserted] = new_mortonbox2discretizationpoints_loc.try_emplace(new_morton_box, values);
            
                    if (!inserted) {
                
                        auto& entry = it->second;

                        entry[0] = std::min(entry[0], values[0]);
                        entry[1] += values[1];

                    }

                }

                mortonbox2discretizationpoints_loc = new_mortonbox2discretizationpoints_loc;

            }   

        }

        void InitializeBoxesAndLevels() 
        {

            std::array<double, 3> min, max;
            ComputeBB(min, max);

            //const double boxsize_x = (max[0] - min[0]) + BBSIZEOFFSET;
            //const double boxsize_y = (max[1] - min[1]) + BBSIZEOFFSET;
            //const double boxsize_z = (max[2] - min[2]) + BBSIZEOFFSET;

            const double boxsize = std::max(max[0] - min[0], std::max(max[1] - min[1], max[2] - min[2])) + BBSIZEOFFSET;
            const double boxsize_x = boxsize;
            const double boxsize_y = boxsize;
            const double boxsize_z = boxsize;

            min[0] -= BBSIZEOFFSET/2.0;
            min[1] -= BBSIZEOFFSET/2.0;
            min[2] -= BBSIZEOFFSET/2.0;

            sorting_.resize(N_loc_orig_);
            std::iota(sorting_.begin(), sorting_.end(), point_low_);  

            if (USE_ADAPTIVITY) {

                ComputeNLevels(min, boxsize_x, boxsize_y, boxsize_z);

            }            

            SortBox(min, boxsize_x / static_cast<double>(1LL << (nlevels_ - 1)), boxsize_y / static_cast<double>(1LL << (nlevels_ - 1)), boxsize_z / static_cast<double>(1LL << (nlevels_ - 1)));

            levels_.resize(nlevels_, nullptr);

            #pragma omp parallel for
            for (int i = 0; i < nlevels_; i++) {

                levels_[i] = new Level(i, min[0], min[1], min[2], boxsize_x / static_cast<double>(1LL << i), boxsize_y / static_cast<double>(1LL << i), boxsize_z / static_cast<double>(1LL << i), wavenumber_);

            }

            InitializeLevelDBoxesData();

            InitializeRelBoxesAllLevels();

            InitializeBoxesDataAllLevels();

        }

        void GetRelevantConeSegmentsDueToCousinSurface(const int level, std::unordered_set<long long>& mortonboxnonrelconesegment)
        {

            const long long ncones = levels_[level]->GetTotalNumberOfCones();

            #pragma omp parallel
            {

            long long old_morton_box = -1;
            std::vector<long long> cousins;
            std::unordered_set<long long> tmpmortonboxnonrelconesegment; 

            #pragma omp for
            for (long long i = 0; i < N_loc_orig_; i++) {

                const double x = x_[i];
                const double y = y_[i];
                const double z = z_[i];

                const long long morton_box = levels_[level]->Point2Morton(x, y, z);

                if (morton_box != old_morton_box) {

                    cousins = levels_[level]->GetCousins(morton_box);
                    old_morton_box = morton_box;

                }

                for (const long long morton_cousin : cousins) {

                    const long long size_cousin = levels_[level]->mortonbox2discretizationpoints_2_[morton_cousin][1];

                    if (USE_ADAPTIVITY && (size_cousin <= MAX_ELEMS_LEAF)) {

                        continue;

                    }

                    const long long nonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_cousin, x, y, z);

                    const long long id = morton_cousin * ncones + nonrel_conesegment;

                    tmpmortonboxnonrelconesegment.insert(id);

                }

            }

            #pragma omp critical
            {

                mortonboxnonrelconesegment.insert(tmpmortonboxnonrelconesegment.begin(), tmpmortonboxnonrelconesegment.end());

            }

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

        void GetRelevantConeSegmentsDueToParents(const int level, std::unordered_set<long long>& mortonboxnonrelconesegment) 
        {

            #pragma omp parallel
            {

            double x, y, z;
            long long old_cocentered_morton_box = -1;
            std::vector<long long> morton_children;
            
            std::unordered_set<long long> mortonboxnonrelconesegment_thread;

            const long long ncones = levels_[level-1]->GetTotalNumberOfCones();
            const long long ncones_2 = levels_[level]->GetTotalNumberOfCones();

            #pragma omp for
            for (const auto & cs : levels_[level-1]->relconesegmentids_) {

                const long long cocentered_morton_box = cs / ncones;
                const long long nonrel_conesegment = cs % ncones;
            
                if (cocentered_morton_box != old_cocentered_morton_box) {

                    morton_children = GetMortonRelChildren(level-1, cocentered_morton_box);

                    if (morton_children.size() == 0) {
                        throw std::logic_error("Every cocentered morton box must have children in get relevant cone segments function.");
                    }

                    old_cocentered_morton_box = cocentered_morton_box;

                }

                for (const auto & morton_child : morton_children) {

                    const std::array<long long, 2> data_child = levels_[level]->mortonbox2discretizationpoints_2_[morton_child];

                    if (USE_ADAPTIVITY && (data_child[1] <= MAX_ELEMS_LEAF)) {

                        continue;

                    }

                    for (int interpiter = 0; interpiter < IPSCHEME_.GetNInterpPoints(); interpiter++) {

                        IPSCHEME_.GetInterpolationPoint(interpiter, x, y, z);

                        levels_[level-1]->Cheb2Cart(cocentered_morton_box, nonrel_conesegment, x, y, z);

                        long long locnonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_child, x, y, z);

                        long long id = morton_child * ncones_2 + locnonrel_conesegment;
                        
                        mortonboxnonrelconesegment_thread.insert(id);

                    } 

                }

            }

            #pragma omp critical
            {
            
                mortonboxnonrelconesegment.insert(mortonboxnonrelconesegment_thread.begin(), mortonboxnonrelconesegment_thread.end());

            }

            mortonboxnonrelconesegment_thread.clear();

            }

        }

        std::vector<long long> MergeRelevantConeSegments(
            const std::unordered_set<long long>& mortonboxnonrelconesegment1,
            const std::unordered_set<long long>& mortonboxnonrelconesegment2) const 
        {

            std::vector<long long> merged_set;
            merged_set.reserve(mortonboxnonrelconesegment1.size() + mortonboxnonrelconesegment2.size());

            merged_set.assign(mortonboxnonrelconesegment1.begin(), mortonboxnonrelconesegment1.end());

            merged_set.insert(merged_set.end(), mortonboxnonrelconesegment2.begin(), mortonboxnonrelconesegment2.end());

            std::sort(merged_set.begin(), merged_set.end());
            auto last_unique = std::unique(merged_set.begin(), merged_set.end());
            merged_set.erase(last_unique, merged_set.end());

            return merged_set;

        }

        void CountingSortDigit_v2(std::vector<long long>& data, long long exp) 
        {

            long long n = data.size();

            std::vector<long long> output(n);
            std::vector<long long> count(10, 0);
            
            for (long long i = 0; i < n; i++) {

                count[(data[i] / exp) % 10]++;

            }
        
            for (int i = 1; i < 10; i++) {

                count[i] += count[i - 1];

            }
        
            for (long long i = n-1; i >= 0; i--) {

                output[count[(data[i] / exp) % 10] - 1] = data[i];
                count[(data[i] / exp) % 10]--;

            }
        
            data = output;

        }

        void LoadBalanceRelevantConeSegments(const int level, std::vector<long long>& mortonboxnonrelconesegment) 
        {

            long long local_max;
            if (mortonboxnonrelconesegment.size() == 0) {
                local_max = -1;
            } else {
                local_max = *std::max_element(mortonboxnonrelconesegment.begin(), mortonboxnonrelconesegment.end());
            }
            long long global_max;
            MPI_Allreduce(&local_max, &global_max, 1, MPI_LONG_LONG, MPI_MAX, mpi_comm_);

            for (long long exp = 1; global_max / exp > 0; exp *= 10) {

                CountingSortDigit_v2(mortonboxnonrelconesegment, exp); 
                auto unique_elems = std::unique(mortonboxnonrelconesegment.begin(), mortonboxnonrelconesegment.end());
                mortonboxnonrelconesegment.erase(unique_elems, mortonboxnonrelconesegment.end());

                std::vector<long long> local_counts(10, 0);  

                for (long long val : mortonboxnonrelconesegment) {

                    local_counts[(val / exp) % 10]++;

                }   
                
                std::vector<long long> global_counts(10, 0);
                MPI_Allreduce(local_counts.data(), global_counts.data(), 10, MPI_LONG_LONG, MPI_SUM, mpi_comm_);

                std::vector<long long> partial_counts(10, 0);
                MPI_Exscan(local_counts.data(), partial_counts.data(), 10, MPI_LONG_LONG, MPI_SUM, mpi_comm_);

                std::vector<long long> global_displs(11, 0);
                for (int i = 1; i < 11; i++) {
                    global_displs[i] = global_displs[i-1] + global_counts[i-1];
                }

                std::vector<MPI_Count> send_counts(comm_size_, 0);
                std::vector<MPI_Aint> send_displs(comm_size_, 0);  
                
                long long N_exp = 0;
                for (int i = 0; i < 10; i++) {
                    N_exp += global_counts[i];
                }

                long long N_loc = N_exp / comm_size_;
                long long N_loc_rem = N_exp % comm_size_;

                std::vector<long long> split_points_loc(comm_size_+1, 0);

                for (int i = 1; i <= comm_size_; i++) {

                    split_points_loc[i] = split_points_loc[i-1] + N_loc;

                    if (N_loc_rem > 0) {

                        split_points_loc[i]++;
                        N_loc_rem--;

                    }

                }

                for (int i = 0; i < comm_size_; ++i) {

                    long long pos_low = split_points_loc[i];
                    long long pos_up = split_points_loc[i+1]-1;
                    
                    auto it_low = std::upper_bound(global_displs.begin(), global_displs.end(), pos_low);
                    int digit_low = static_cast<int>(std::distance(global_displs.begin(), it_low)) - 1;
                    if (digit_low < 0) { 
                        digit_low = 0;            
                    }
                    if (digit_low >= 10) {                     
                        digit_low = 9; 
                    }

                    auto it_up = std::upper_bound(global_displs.begin(), global_displs.end(), pos_up);
                    int digit_up = static_cast<int>(std::distance(global_displs.begin(), it_up)) - 1;
                    if (digit_up < 0) { 
                        digit_up = 0;            
                    }
                    if (digit_up >= 10) {                     
                        digit_up = 9; 
                    }

                    for (int j = digit_low; j <= digit_up; j++) {

                        if (local_counts[j] == 0) {
                            continue;
                        }

                        long long needed_count;

                        if (digit_low == digit_up) {
                            needed_count = pos_up + 1 - pos_low;
                        } else if (j == digit_low) {
                            needed_count = global_displs[j+1] - pos_low;
                        } else if (j == digit_up) {
                            needed_count = pos_up + 1 - global_displs[j];
                        } else {
                            needed_count = global_displs[j+1] - global_displs[j];
                        }

                        if (partial_counts[j] >= needed_count) {
                            partial_counts[j] -= needed_count;
                        } else {
                            long long total_present = local_counts[j];
                            long long already_sent = partial_counts[j];
                            long long need_to_be_sent = needed_count - already_sent;
                            long long total_sent = std::min(total_present, need_to_be_sent); 
                            send_counts[i] += total_sent;
                            local_counts[j] -= total_sent; 
                            partial_counts[j] = 0;
                        }                           

                    }

                    if (i > 0) {

                        send_displs[i] = send_displs[i-1] + send_counts[i-1];

                    }

                }

                std::vector<MPI_Count> recv_counts(comm_size_, 0);
                std::vector<MPI_Aint> recv_displs(comm_size_, 0);

                MPI_Alltoall(send_counts.data(), 1, MPI_COUNT, recv_counts.data(), 1, MPI_COUNT, mpi_comm_);
                
                for (int i = 1; i < comm_size_; ++i) {

                    recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];

                }

                std::vector<long long> new_mortonboxnonrelconesegment(split_points_loc[comm_rank_+1] - split_points_loc[comm_rank_]);

                MPI_Alltoallv_c(mortonboxnonrelconesegment.data(), send_counts.data(), send_displs.data(), MPI_LONG_LONG,
                                new_mortonboxnonrelconesegment.data(), recv_counts.data(), recv_displs.data(), MPI_LONG_LONG,
                                mpi_comm_);

                mortonboxnonrelconesegment = new_mortonboxnonrelconesegment;

                CountingSortDigit_v2(mortonboxnonrelconesegment, exp); 
                auto unique_elems_2 = std::unique(mortonboxnonrelconesegment.begin(), mortonboxnonrelconesegment.end());
                mortonboxnonrelconesegment.erase(unique_elems_2, mortonboxnonrelconesegment.end());

            }  

            long long first_elem;
            long long last_elem;

            if (mortonboxnonrelconesegment.size() == 0) {
                first_elem = -1;
                last_elem = -1;
            } else {
                first_elem = mortonboxnonrelconesegment[0];
                last_elem = mortonboxnonrelconesegment.back();
            }

            std::vector<long long> first_elem_all(comm_size_);
            std::vector<long long> last_elem_all(comm_size_);

            MPI_Allgather(&first_elem, 1, MPI_LONG_LONG, &first_elem_all[0], 1, MPI_LONG_LONG, mpi_comm_);
            MPI_Allgather(&last_elem, 1, MPI_LONG_LONG, &last_elem_all[0], 1, MPI_LONG_LONG, mpi_comm_);

            first_elem_all.push_back(last_elem_all.back()+1);

            if (first_elem_all.back() == 0) {

                first_elem_all[comm_size_] = -1;

                auto it = std::find(first_elem_all.begin(), first_elem_all.end(), -1);
                int pos = std::distance(first_elem_all.begin(), it);

                for (int rank = pos; rank <= comm_size_; rank++) {

                    first_elem_all[rank] = first_elem_all[pos-1] + 1;

                }

            }

            levels_[level]->split_points_relconesegments_ = std::move(first_elem_all);
            levels_[level]->relconesegmentids_ = std::move(mortonboxnonrelconesegment);

        }

        void SetRelevantConeSegmentsNotInRank(int level, 
                                              const std::unordered_set<long long>& interpolation_conesegments,
                                              const std::unordered_set<long long>& propagation_conesegments) 
        {

            std::vector<std::vector<long long>> propagation_not_in_rank(comm_size_);

            for (long long cone_segment : propagation_conesegments) {

                auto it = std::upper_bound(levels_[level]->split_points_relconesegments_.begin(), levels_[level]->split_points_relconesegments_.end(), cone_segment);

                int rank = static_cast<int>(std::distance(levels_[level]->split_points_relconesegments_.begin(), it)) - 1;
                if (rank < 0) rank = 0;
                if (rank >= comm_size_) rank = comm_size_ - 1;            

                if (rank != comm_rank_) {

                    propagation_not_in_rank[rank].push_back(cone_segment);

                }

            }

            levels_[level]->send_counts_propagation_.resize(comm_size_);
            levels_[level]->sdispls_propagation_.resize(comm_size_);
            MPI_Count total_send_local = 0;

            for (int i = 0; i < comm_size_; ++i) {

                levels_[level]->sdispls_propagation_[i] = total_send_local;
                levels_[level]->send_counts_propagation_[i] = static_cast<MPI_Count>(propagation_not_in_rank[i].size());
                total_send_local += levels_[level]->send_counts_propagation_[i];

            }

            levels_[level]->flat_send_buffer_propagation_.resize(static_cast<size_t>(total_send_local));

            MPI_Aint current_pos = 0;
    
            for (int i = 0; i < comm_size_; ++i) {

                std::copy(propagation_not_in_rank[i].begin(), propagation_not_in_rank[i].end(),
                          levels_[level]->flat_send_buffer_propagation_.begin() + current_pos);
                current_pos += levels_[level]->send_counts_propagation_[i];

            }

            propagation_not_in_rank.clear();

            std::vector<MPI_Count> recv_counts_propagation(comm_size_); 

            MPI_Alltoall(levels_[level]->send_counts_propagation_.data(), 1, MPI_COUNT, recv_counts_propagation.data(), 1, MPI_COUNT, mpi_comm_);

            std::vector<MPI_Aint> rdispls_propagation(comm_size_);
            MPI_Count total_received_size = 0;

            for (int i = 0; i < comm_size_; ++i) {

                rdispls_propagation[i] = total_received_size;
                total_received_size += recv_counts_propagation[i];

            }

            std::vector<long long> flat_recv_buffer_propagation(static_cast<size_t>(total_received_size));

            MPI_Alltoallv_c(levels_[level]->flat_send_buffer_propagation_.data(), levels_[level]->send_counts_propagation_.data(), levels_[level]->sdispls_propagation_.data(), MPI_LONG_LONG,
                            flat_recv_buffer_propagation.data(), recv_counts_propagation.data(), rdispls_propagation.data(), MPI_LONG_LONG, 
                            mpi_comm_);   

            std::vector<long long> pos_send_back_propagation(total_received_size);

            #pragma omp parallel for
            for (long long i = 0; i < total_received_size; i++) {
            
                const long long cs = flat_recv_buffer_propagation[i];
                
                auto it = std::lower_bound(levels_[level]->relconesegmentids_.begin(), levels_[level]->relconesegmentids_.end(), cs);
                
                const long long idx = std::distance(levels_[level]->relconesegmentids_.begin(), it);
                
                pos_send_back_propagation[i] = idx * P_;
            
            }
            
            flat_recv_buffer_propagation.clear();
            
            std::vector<long long> pos_recv_back_propagation(total_send_local);
            
            MPI_Alltoallv_c(pos_send_back_propagation.data(), recv_counts_propagation.data(), rdispls_propagation.data(), MPI_LONG_LONG,
                            pos_recv_back_propagation.data(), levels_[level]->send_counts_propagation_.data(), levels_[level]->sdispls_propagation_.data(), MPI_LONG_LONG, 
                            mpi_comm_); 

            pos_send_back_propagation.clear();
            
            levels_[level]->flat_pos_propagation_ = std::move(pos_recv_back_propagation);

            std::vector<std::vector<long long>> interpolation_not_in_rank(comm_size_);

            for (long long cone_segment : interpolation_conesegments) {

                auto it = std::upper_bound(levels_[level]->split_points_relconesegments_.begin(), levels_[level]->split_points_relconesegments_.end(), cone_segment);

                int rank = static_cast<int>(std::distance(levels_[level]->split_points_relconesegments_.begin(), it)) - 1;
                if (rank < 0) rank = 0;
                if (rank >= comm_size_) rank = comm_size_ - 1;            

                if (rank != comm_rank_) {

                    interpolation_not_in_rank[rank].push_back(cone_segment);

                }

            }

            levels_[level]->send_counts_interpolation_.resize(comm_size_);
            levels_[level]->sdispls_interpolation_.resize(comm_size_);
            total_send_local = 0;

            for (int i = 0; i < comm_size_; ++i) {

                levels_[level]->sdispls_interpolation_[i] = total_send_local;
                levels_[level]->send_counts_interpolation_[i] = static_cast<MPI_Count>(interpolation_not_in_rank[i].size());
                total_send_local += levels_[level]->send_counts_interpolation_[i];

            }

            levels_[level]->flat_send_buffer_interpolation_.resize(static_cast<size_t>(total_send_local));

            current_pos = 0;
    
            for (int i = 0; i < comm_size_; ++i) {

                std::copy(interpolation_not_in_rank[i].begin(), interpolation_not_in_rank[i].end(),
                          levels_[level]->flat_send_buffer_interpolation_.begin() + current_pos);
                current_pos += levels_[level]->send_counts_interpolation_[i];

            }

            interpolation_not_in_rank.clear();

            std::vector<MPI_Count> recv_counts_interpolation(comm_size_); 

            MPI_Alltoall(levels_[level]->send_counts_interpolation_.data(), 1, MPI_COUNT, recv_counts_interpolation.data(), 1, MPI_COUNT, mpi_comm_);

            std::vector<MPI_Aint> rdispls_interpolation(comm_size_);
            total_received_size = 0;

            for (int i = 0; i < comm_size_; ++i) {

                rdispls_interpolation[i] = total_received_size;
                total_received_size += recv_counts_interpolation[i];

            }

            std::vector<long long> flat_recv_buffer_interpolation(static_cast<size_t>(total_received_size));

            MPI_Alltoallv_c(levels_[level]->flat_send_buffer_interpolation_.data(), levels_[level]->send_counts_interpolation_.data(), levels_[level]->sdispls_interpolation_.data(), MPI_LONG_LONG,
                            flat_recv_buffer_interpolation.data(), recv_counts_interpolation.data(), rdispls_interpolation.data(), MPI_LONG_LONG, 
                            mpi_comm_);   

            std::vector<long long> pos_send_back_interpolation(total_received_size);

            #pragma omp parallel for
            for (long long i = 0; i < total_received_size; i++) {
            
                const long long cs = flat_recv_buffer_interpolation[i];
                
                auto it = std::lower_bound(levels_[level]->relconesegmentids_.begin(), levels_[level]->relconesegmentids_.end(), cs);
                
                const long long idx = std::distance(levels_[level]->relconesegmentids_.begin(), it);
                
                pos_send_back_interpolation[i] = idx * P_;
            
            }
            
            flat_recv_buffer_interpolation.clear();
            
            std::vector<long long> pos_recv_back_interpolation(total_send_local);
            
            MPI_Alltoallv_c(pos_send_back_interpolation.data(), recv_counts_interpolation.data(), rdispls_interpolation.data(), MPI_LONG_LONG,
                            pos_recv_back_interpolation.data(), levels_[level]->send_counts_interpolation_.data(), levels_[level]->sdispls_interpolation_.data(), MPI_LONG_LONG, 
                            mpi_comm_); 

            pos_send_back_interpolation.clear();
            
            levels_[level]->flat_pos_interpolation_ = std::move(pos_recv_back_interpolation);

        }

        void SetRelevantChildrenOfCocenteredBoxes(int level)
        {

            if (level >= nlevels_ - 1) { 
                
                return;
            
            }

            std::vector<long long> morton_box_loc;
            std::vector<long long> position_loc;
            std::vector<long long> size_loc;

            morton_box_loc.reserve(levels_[level+1]->mortonbox2discretizationpoints_.size());
            position_loc.reserve(levels_[level+1]->mortonbox2discretizationpoints_.size());
            size_loc.reserve(levels_[level+1]->mortonbox2discretizationpoints_.size());

            for (const auto & [morton_box, values] : levels_[level+1]->mortonbox2discretizationpoints_) {

                morton_box_loc.push_back(morton_box);
                position_loc.push_back(values[0]);
                size_loc.push_back(values[1]);

            }

            MPI_Count total_size_loc = morton_box_loc.size();

            std::vector<MPI_Count> recv_counts(comm_size_);
            MPI_Allgather(&total_size_loc, 1, MPI_COUNT, &recv_counts[0], 1, MPI_COUNT, mpi_comm_);

            std::vector<MPI_Aint> displs(comm_size_);
            MPI_Count total_size = 0;

            for (int i = 0; i < comm_size_; i++) {

                displs[i] = total_size;
                total_size += recv_counts[i];

            }

            std::vector<long long> morton_box_all(static_cast<size_t>(total_size));
            std::vector<long long> position_all(static_cast<size_t>(total_size));
            std::vector<long long> size_all(static_cast<size_t>(total_size));

            MPI_Allgatherv_c(&morton_box_loc[0], total_size_loc, MPI_LONG_LONG, &morton_box_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
            MPI_Allgatherv_c(&position_loc[0], total_size_loc, MPI_LONG_LONG, &position_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
            MPI_Allgatherv_c(&size_loc[0], total_size_loc, MPI_LONG_LONG, &size_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);

            std::vector<long long>().swap(morton_box_loc);
            std::vector<long long>().swap(position_loc);
            std::vector<long long>().swap(size_loc);

            std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_all;

            for (long long i = 0; i < total_size; i++) {

                const long long morton_box = morton_box_all[i];
                const std::array<long long, 2> value{position_all[i], size_all[i]};

                auto [it, inserted] = mortonbox2discretizationpoints_all.try_emplace(morton_box, value);

                if (!inserted) {

                    it->second[0] = std::min(it->second[0], value[0]);
                    it->second[1] += value[1];

                }

            }                

            std::vector<long long>().swap(morton_box_all);
            std::vector<long long>().swap(position_all);
            std::vector<long long>().swap(size_all);

            const long long ncones = levels_[level]->GetTotalNumberOfCones();

            for (const auto & idx : levels_[level]->relconesegmentids_) {

                long long morton_box = idx / ncones;

                for (long long childiter = morton_box*8; childiter < (morton_box+1)*8; childiter++) {

                    if (mortonbox2discretizationpoints_all.count(childiter) != 0 && levels_[level+1]->mortonidofrelboxes_2_.count(childiter) == 0) {
                            
                        levels_[level+1]->mortonidofrelboxes_2_.insert(childiter);
                        levels_[level+1]->mortonbox2discretizationpoints_2_[childiter] = mortonbox2discretizationpoints_all[childiter];

                    }

                }

            }

            mortonbox2discretizationpoints_all.clear();

        }      

        void UpdateBoxesData(int level) 
        {

            std::vector<long long> morton_box_loc;
            std::vector<long long> position_loc;
            std::vector<long long> size_loc;

            morton_box_loc.reserve(levels_[level]->mortonbox2discretizationpoints_.size());
            position_loc.reserve(levels_[level]->mortonbox2discretizationpoints_.size());
            size_loc.reserve(levels_[level]->mortonbox2discretizationpoints_.size());

            for (const auto & [morton_box, values] : levels_[level]->mortonbox2discretizationpoints_) {

                morton_box_loc.push_back(morton_box);
                position_loc.push_back(values[0]);
                size_loc.push_back(values[1]);

            }

            MPI_Count total_size_loc = morton_box_loc.size();

            std::vector<MPI_Count> recv_counts(comm_size_);
            MPI_Allgather(&total_size_loc, 1, MPI_COUNT, &recv_counts[0], 1, MPI_COUNT, mpi_comm_);

            std::vector<MPI_Aint> displs(comm_size_);
            MPI_Count total_size = 0;

            for (int i = 0; i < comm_size_; i++) {

                displs[i] = total_size;
                total_size += recv_counts[i];

            }

            std::vector<long long> morton_box_all(static_cast<size_t>(total_size));
            std::vector<long long> position_all(static_cast<size_t>(total_size));
            std::vector<long long> size_all(static_cast<size_t>(total_size));

            MPI_Allgatherv_c(&morton_box_loc[0], total_size_loc, MPI_LONG_LONG, &morton_box_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
            MPI_Allgatherv_c(&position_loc[0], total_size_loc, MPI_LONG_LONG, &position_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
            MPI_Allgatherv_c(&size_loc[0], total_size_loc, MPI_LONG_LONG, &size_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);

            std::vector<long long>().swap(morton_box_loc);
            std::vector<long long>().swap(position_loc);
            std::vector<long long>().swap(size_loc);
            
            std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_all;

            for (long long i = 0; i < total_size; i++) {

                const long long morton_box = morton_box_all[i];
                const std::array<long long, 2> value{position_all[i], size_all[i]};

                auto [it, inserted] = mortonbox2discretizationpoints_all.try_emplace(morton_box, value);

                if (!inserted) {

                    it->second[0] = std::min(it->second[0], value[0]);
                    it->second[1] += value[1];

                }

            } 

            std::vector<long long>().swap(morton_box_all);
            std::vector<long long>().swap(position_all);
            std::vector<long long>().swap(size_all);

            const long long ncones = levels_[level]->GetTotalNumberOfCones();

            for (const auto & idx : levels_[level]->relconesegmentids_) {

                long long morton_box = idx / ncones;
                
                if (levels_[level]->mortonbox2discretizationpoints_2_.count(morton_box) == 0) {

                    levels_[level]->mortonidofrelboxes_2_.insert(morton_box);
                    levels_[level]->mortonbox2discretizationpoints_2_[morton_box] = mortonbox2discretizationpoints_all[morton_box];
                
                }

            }

            mortonbox2discretizationpoints_all.clear();

        }

        void InitializeRelevantConeSegments() 
        {

            std::unordered_set<long long> interpolation_conesegments;
            std::unordered_set<long long> propagation_conesegments;
 
            for (int level = 2; level < nlevels_; level++) {    

                GetRelevantConeSegmentsDueToCousinSurface(level, interpolation_conesegments);

                if (level > 2) {

                    GetRelevantConeSegmentsDueToParents(level, propagation_conesegments);

                }

                std::vector<long long> all_conesegments = MergeRelevantConeSegments(interpolation_conesegments, propagation_conesegments);

                LoadBalanceRelevantConeSegments(level, all_conesegments);
                
                std::vector<long long>().swap(all_conesegments);

                SetRelevantConeSegmentsNotInRank(level, interpolation_conesegments, propagation_conesegments);

                interpolation_conesegments.clear();
                propagation_conesegments.clear();

                SetRelevantChildrenOfCocenteredBoxes(level);

                UpdateBoxesData(level);

            } 

        }  

        void InitializePointsLevelDNeighbours(std::vector<std::unordered_set<long long>>& points_not_in_rank)
        {

            std::vector<long long> vec_mortonids(levels_.back()->mortonidofrelboxes_.begin(), levels_.back()->mortonidofrelboxes_.end());

            #pragma omp parallel
            {

            std::vector<std::unordered_set<long long>> points_not_in_rank_thread(comm_size_);
            
            #pragma omp for
            for (long long i = 0; i < vec_mortonids.size(); i++) {

                const long long morton_box = vec_mortonids[i];

                std::vector<long long> neighbours = levels_.back()->GetNeighbours(morton_box);

                for (const auto & neighbour_morton : neighbours) {

                    const std::array<long long, 2> points_data = levels_.back()->mortonbox2discretizationpoints_2_[neighbour_morton];

                    const long long points_begin = points_data[0];
                    const long long npoints = points_data[1];

                    for (long long j = 0; j < npoints; j++) {

                        auto it = std::upper_bound(split_points_orig_.begin(), split_points_orig_.end(), points_begin + j);
                        int rank = std::clamp<int>(static_cast<int>(std::distance(split_points_orig_.begin(), it)) - 1, 0, comm_size_-1);

                        if (rank != comm_rank_) points_not_in_rank_thread[rank].insert(points_begin + j);

                    }

                }

            }

            #pragma omp critical
            {
            
                for (int i = 0; i < comm_size_; i++) {

                    points_not_in_rank[i].insert(points_not_in_rank_thread[i].begin(), points_not_in_rank_thread[i].end());

                }
                
            }

            }

        }

        void InitializePointsLevelDRelCS(std::vector<std::unordered_set<long long>>& points_not_in_rank)
        {

            std::vector<long long> vec_mortonids(levels_.back()->relconesegmentids_.begin(), levels_.back()->relconesegmentids_.end());
            const long long ncones = levels_.back()->GetTotalNumberOfCones();
            
            #pragma omp parallel
            {

            std::vector<std::unordered_set<long long>> points_not_in_rank_thread(comm_size_);
            long long old_morton_box = -1;
            
            #pragma omp for
            for (long long i = 0; i < vec_mortonids.size(); i++) {
                
                const long long morton_box = vec_mortonids[i] / ncones;

                if (morton_box == old_morton_box) continue;

                old_morton_box = morton_box;

                const std::array<long long, 2> points_data = levels_.back()->mortonbox2discretizationpoints_2_[morton_box];

                const long long points_begin = points_data[0];
                const long long npoints = points_data[1];

                for (long long j = 0; j < npoints; j++) {

                    auto it = std::upper_bound(split_points_orig_.begin(), split_points_orig_.end(), points_begin + j);
                    int rank = std::clamp<int>(static_cast<int>(std::distance(split_points_orig_.begin(), it)) - 1, 0, comm_size_-1);

                    if (rank != comm_rank_) points_not_in_rank_thread[rank].insert(points_begin + j);

                }

            }

            #pragma omp critical
            {
            
                for (int i = 0; i < comm_size_; i++) {

                    points_not_in_rank[i].insert(points_not_in_rank_thread[i].begin(), points_not_in_rank_thread[i].end());

                }
                
            }

            }

        }

        void InitializePointsLeveldCousins(int level,
                                           std::vector<std::unordered_set<long long>>& points_not_in_rank)
        {

            std::vector<long long> vec_mortonids(levels_[level]->mortonidofrelboxes_.begin(), levels_[level]->mortonidofrelboxes_.end());

            #pragma omp parallel
            {

            std::vector<std::unordered_set<long long>> points_not_in_rank_thread(comm_size_);
            
            #pragma omp for
            for (long long i = 0; i < vec_mortonids.size(); i++) {

                const long long morton_box = vec_mortonids[i];

                std::vector<long long> cousins = levels_[level]->GetCousins(morton_box);

                for (const auto & cousin_morton : cousins) {

                    const std::array<long long, 2> points_data = levels_[level]->mortonbox2discretizationpoints_2_[cousin_morton];

                    long long points_begin = points_data[0];
                    long long npoints = points_data[1];

                    if (USE_ADAPTIVITY && npoints <= MAX_ELEMS_LEAF) {

                        for (long long j = 0; j < npoints; j++) {

                            auto it = std::upper_bound(split_points_orig_.begin(), split_points_orig_.end(), points_begin + j);
                            int rank = std::clamp<int>(static_cast<int>(std::distance(split_points_orig_.begin(), it)) - 1, 0, comm_size_-1);
    
                            if (rank != comm_rank_) points_not_in_rank_thread[rank].insert(points_begin + j);
    
                        }

                    }

                }

            }

            #pragma omp critical
            {
            
                for (int i = 0; i < comm_size_; i++) {

                    points_not_in_rank[i].insert(points_not_in_rank_thread[i].begin(), points_not_in_rank_thread[i].end());

                }
                
            }

            }

        }

        void InitializePointsLeveldRelCS(int level,
                                         std::vector<std::unordered_set<long long>>& points_not_in_rank)
        {

            std::vector<long long> vec_mortonids(levels_[level-1]->relconesegmentids_.begin(), levels_[level-1]->relconesegmentids_.end());
            const long long ncones = levels_[level-1]->GetTotalNumberOfCones();

            #pragma omp parallel
            {

            std::vector<std::unordered_set<long long>> points_not_in_rank_thread(comm_size_);
            long long old_morton_box = -1;
            
            #pragma omp for
            for (long long i = 0; i < vec_mortonids.size(); i++) {

                long long morton_box = vec_mortonids[i] / ncones;

                if (morton_box == old_morton_box) continue;

                old_morton_box = morton_box;

                std::vector<long long> morton_children = GetMortonRelChildren(level-1, morton_box);

                for (const long long morton_child : morton_children) {

                    const std::array<long long, 2> points_data = levels_[level]->mortonbox2discretizationpoints_2_[morton_child];

                    long long points_begin = points_data[0];
                    long long npoints = points_data[1];

                    if (USE_ADAPTIVITY && npoints <= MAX_ELEMS_LEAF) {

                        for (long long j = 0; j < npoints; j++) {

                            auto it = std::upper_bound(split_points_orig_.begin(), split_points_orig_.end(), points_begin + j);
                            int rank = std::clamp<int>(static_cast<int>(std::distance(split_points_orig_.begin(), it)) - 1, 0, comm_size_-1);
    
                            if (rank != comm_rank_) points_not_in_rank_thread[rank].insert(points_begin + j);
    
                        }

                    }

                }

            }

            #pragma omp critical
            {
            
                for (int i = 0; i < comm_size_; i++) {

                    points_not_in_rank[i].insert(points_not_in_rank_thread[i].begin(), points_not_in_rank_thread[i].end());

                }
                
            }

            }

        }

        void GetDataNotInRank(std::vector<std::unordered_set<long long>>& points_not_in_rank)
        {

            std::vector<MPI_Count> send_counts_points(comm_size_);
            std::vector<MPI_Aint> sdispls_points(comm_size_);
            MPI_Count total_send_points = 0;    

            for (int i = 0; i < comm_size_; ++i) {

                send_counts_points[i] = static_cast<MPI_Count>(points_not_in_rank[i].size());
                sdispls_points[i] = total_send_points;
                total_send_points += send_counts_points[i];

            }

            std::vector<long long> flat_send_points_buffer(static_cast<size_t>(total_send_points));
            
            MPI_Aint current_point_pos = 0;
    
            for (int i = 0; i < comm_size_; ++i) {

                std::copy(points_not_in_rank[i].begin(), points_not_in_rank[i].end(),
                          flat_send_points_buffer.begin() + current_point_pos);
                current_point_pos += send_counts_points[i];

            }

            std::vector<MPI_Count> recv_counts_points(comm_size_);
    
            MPI_Alltoall(send_counts_points.data(), 1, MPI_COUNT, recv_counts_points.data(), 1, MPI_COUNT, mpi_comm_);

            std::vector<MPI_Aint> rdispls_points(comm_size_);
            MPI_Count total_recv_points = 0;

            for (int i = 0; i < comm_size_; ++i) {

                rdispls_points[i] = total_recv_points;
                total_recv_points += recv_counts_points[i];

            }

            std::vector<long long> flat_recv_points_buffer(static_cast<size_t>(total_recv_points));
            
            MPI_Alltoallv_c(flat_send_points_buffer.data(), send_counts_points.data(), sdispls_points.data(), MPI_LONG_LONG,
                            flat_recv_points_buffer.data(), recv_counts_points.data(), rdispls_points.data(), MPI_LONG_LONG,
                            mpi_comm_);

            std::vector<MPI_Count> send_counts_points_back(comm_size_);
            std::vector<MPI_Aint> sdispls_points_back(comm_size_);
            MPI_Count total_send_points_back = 0;

            for (int i = 0; i < comm_size_; ++i) {

                send_counts_points_back[i] = recv_counts_points[i] * 7;
                sdispls_points_back[i] = total_send_points_back;
                total_send_points_back += send_counts_points_back[i];

            }

            std::vector<double> flat_send_points_back_buffer(static_cast<size_t>(total_send_points_back));

            for (int sender_rank = 0; sender_rank < comm_size_; ++sender_rank) {
        
                MPI_Aint recv_pos = rdispls_points[sender_rank];
                MPI_Aint out_pos = sdispls_points_back[sender_rank];

                const MPI_Count num_points = recv_counts_points[sender_rank];

                for (MPI_Count i = 0; i < num_points; ++i) {
            
                    const long long point = flat_recv_points_buffer[recv_pos + i] - point_low_;

                    flat_send_points_back_buffer[out_pos + 7*i + 0] = x_[point];
                    flat_send_points_back_buffer[out_pos + 7*i + 1] = y_[point];
                    flat_send_points_back_buffer[out_pos + 7*i + 2] = z_[point];
                    flat_send_points_back_buffer[out_pos + 7*i + 3] = normal_x_[point];
                    flat_send_points_back_buffer[out_pos + 7*i + 4] = normal_y_[point];
                    flat_send_points_back_buffer[out_pos + 7*i + 5] = normal_z_[point];
                    flat_send_points_back_buffer[out_pos + 7*i + 6] = sorting_[point];

                }

            }            

            std::vector<long long>().swap(flat_recv_points_buffer);

            std::vector<MPI_Count> recv_counts_points_back(comm_size_);
    
            MPI_Alltoall(send_counts_points_back.data(), 1, MPI_COUNT,
                         recv_counts_points_back.data(), 1, MPI_COUNT, mpi_comm_);

            MPI_Count total_recv_points_back = 0;
            std::vector<MPI_Aint> rdispls_points_back(comm_size_);

            for (int i = 0; i < comm_size_; ++i) {

                rdispls_points_back[i] = total_recv_points_back;
                total_recv_points_back += recv_counts_points_back[i];

            }

            std::vector<double> flat_recv_points_buffer_back(static_cast<size_t>(total_recv_points_back));

            MPI_Alltoallv_c(flat_send_points_back_buffer.data(), send_counts_points_back.data(), sdispls_points_back.data(), MPI_DOUBLE,
                            flat_recv_points_buffer_back.data(), recv_counts_points_back.data(), rdispls_points_back.data(), MPI_DOUBLE,
                            mpi_comm_);
    
            std::vector<double>().swap(flat_send_points_back_buffer);

            for (MPI_Count i = 0; i < total_send_points; i++) {

                const long long point = flat_send_points_buffer[i];

                const double x = flat_recv_points_buffer_back[7 * i];
                const double y = flat_recv_points_buffer_back[7 * i + 1];
                const double z = flat_recv_points_buffer_back[7 * i + 2];

                const double nx = flat_recv_points_buffer_back[7 * i + 3];
                const double ny = flat_recv_points_buffer_back[7 * i + 4];
                const double nz = flat_recv_points_buffer_back[7 * i + 5];

                const double sorting = flat_recv_points_buffer_back[7 * i + 6];

                data_not_in_rank_[point] = {x, y, z, nx, ny, nz, sorting};

            }

            std::vector<long long>().swap(flat_send_points_buffer);
            std::vector<double>().swap(flat_recv_points_buffer_back);

        }

        void InitializePointsNotInRank() 
        {

            std::vector<std::unordered_set<long long>> points_not_in_rank(comm_size_);

            InitializePointsLevelDNeighbours(points_not_in_rank);

            InitializePointsLevelDRelCS(points_not_in_rank);

            for (int level = nlevels_-1; level >= 2; level--) {

                if (level > 2) {

                    InitializePointsLeveldRelCS(level, points_not_in_rank);

                }

                InitializePointsLeveldCousins(level, points_not_in_rank);

            }

            GetDataNotInRank(points_not_in_rank);

            points_not_in_rank_ = std::move(points_not_in_rank);

        }

        void GetDensityNotInRank(const std::vector<std::complex<double>>& density_loc,
                                 std::unordered_map<long long, std::complex<double>>& density_not_in_rank)
        {           

            std::vector<MPI_Count> send_counts_points(comm_size_);
            std::vector<MPI_Aint> sdispls_points(comm_size_);
            MPI_Count total_send_points = 0;    

            for (int i = 0; i < comm_size_; ++i) {

                send_counts_points[i] = static_cast<MPI_Count>(points_not_in_rank_[i].size());
                sdispls_points[i] = total_send_points;
                total_send_points += send_counts_points[i];

            }

            density_not_in_rank.reserve(static_cast<size_t>(total_send_points));

            std::vector<long long> flat_send_points_buffer(static_cast<size_t>(total_send_points));
            
            MPI_Aint current_point_pos = 0;
    
            for (int i = 0; i < comm_size_; ++i) {

                std::copy(points_not_in_rank_[i].begin(), points_not_in_rank_[i].end(),
                          flat_send_points_buffer.begin() + current_point_pos);
                current_point_pos += send_counts_points[i];

            }

            std::vector<MPI_Count> recv_counts_points(comm_size_);
    
            MPI_Alltoall(send_counts_points.data(), 1, MPI_COUNT, recv_counts_points.data(), 1, MPI_COUNT, mpi_comm_);

            std::vector<MPI_Aint> rdispls_points(comm_size_);
            MPI_Count total_recv_points = 0;

            for (int i = 0; i < comm_size_; ++i) {

                rdispls_points[i] = total_recv_points;
                total_recv_points += recv_counts_points[i];

            }

            std::vector<long long> flat_recv_points_buffer(static_cast<size_t>(total_recv_points));
            
            MPI_Alltoallv_c(flat_send_points_buffer.data(), send_counts_points.data(), sdispls_points.data(), MPI_LONG_LONG,
                            flat_recv_points_buffer.data(), recv_counts_points.data(), rdispls_points.data(), MPI_LONG_LONG,
                            mpi_comm_);

            std::vector<MPI_Count> send_counts_points_back(comm_size_);
            std::vector<MPI_Aint> sdispls_points_back(comm_size_);
            MPI_Count total_send_points_back = 0;

            for (int i = 0; i < comm_size_; ++i) {

                send_counts_points_back[i] = recv_counts_points[i];
                sdispls_points_back[i] = total_send_points_back;
                total_send_points_back += send_counts_points_back[i];

            }

            std::vector<std::complex<double>> flat_send_points_back_buffer(static_cast<size_t>(total_send_points_back));

            for (int sender_rank = 0; sender_rank < comm_size_; ++sender_rank) {
        
                MPI_Aint recv_pos = rdispls_points[sender_rank];
                MPI_Aint out_pos = sdispls_points_back[sender_rank];

                const MPI_Count num_points = recv_counts_points[sender_rank];

                for (MPI_Count i = 0; i < num_points; ++i) {
            
                    const long long point = flat_recv_points_buffer[recv_pos + i] - point_low_;

                    flat_send_points_back_buffer[out_pos + i] = density_loc[point];

                }

            }            

            std::vector<long long>().swap(flat_recv_points_buffer); 

            std::vector<MPI_Count> recv_counts_points_back(comm_size_);
    
            MPI_Alltoall(send_counts_points_back.data(), 1, MPI_COUNT,
                         recv_counts_points_back.data(), 1, MPI_COUNT, mpi_comm_);

            MPI_Count total_recv_points_back = 0;
            std::vector<MPI_Aint> rdispls_points_back(comm_size_);

            for (int i = 0; i < comm_size_; ++i) {

                rdispls_points_back[i] = total_recv_points_back;
                total_recv_points_back += recv_counts_points_back[i];

            }

            std::vector<std::complex<double>> flat_recv_points_buffer_back(static_cast<size_t>(total_recv_points_back));

            MPI_Alltoallv_c(flat_send_points_back_buffer.data(), send_counts_points_back.data(), sdispls_points_back.data(), MPI_DOUBLE_COMPLEX,
                            flat_recv_points_buffer_back.data(), recv_counts_points_back.data(), rdispls_points_back.data(), MPI_DOUBLE_COMPLEX,
                            mpi_comm_);
    
            std::vector<std::complex<double>>().swap(flat_send_points_back_buffer);

            for (MPI_Count i = 0; i < total_send_points; i++) {

                const long long point = flat_send_points_buffer[i];

                const std::complex<double> value = flat_recv_points_buffer_back[i];

                density_not_in_rank[point] = value;

            }

            std::vector<long long>().swap(flat_send_points_buffer);
            std::vector<std::complex<double>>().swap(flat_recv_points_buffer_back);

        }

        void ZeroSolution(const long long N) 
        {

            solution_.assign(N, {0.0, 0.0});

        }

        template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
                               const std::complex<double>, const std::complex<double>, const std::complex<double>, 
                               std::complex<double>&)>
        void SingularInteractions(const std::vector<std::complex<double>>& density,
                                  const std::unordered_map<long long, std::complex<double>>& density_not_in_rank,
                                  const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
                                  const bool exclude_pts_RP)
        {

            const long long N = x.size();
            
            const long long patch_size = N_PTS_PER_PATCH[0] * N_PTS_PER_PATCH[1];

            #pragma omp parallel
            {

            long long old_mortonbox = -1;
            std::vector<long long> neighbours;

            #pragma omp for schedule(dynamic)
            for (long long i = 0; i < N; i++) {

                long long point;
                if (exclude_pts_RP) point = sorting_[i];

                const double xp = x[i];
                const double yp = y[i];
                const double zp = z[i];

                const long long mortonbox = levels_.back()->Point2Morton(xp, yp, zp); 

                if (mortonbox != old_mortonbox) {

                    neighbours = levels_.back()->GetNeighbours(mortonbox);

                    old_mortonbox = mortonbox;

                }

                for (const long long neighbourmorton : neighbours) {

                    const std::array<long long, 2> points_data_neig = levels_.back()->mortonbox2discretizationpoints_2_.at(neighbourmorton);

                    long long points_begin = points_data_neig[0];
                    long long npoints = points_data_neig[1];

                    const long long points_end = points_begin + npoints;

                    long long k_local_start = 0, k_local_end = 0;

                    const long long intersection_start = std::max(points_begin, point_low_);
                    const long long intersection_end = std::min(points_end, point_up_);

                    if (intersection_start < intersection_end) {
                        
                        k_local_start = intersection_start - points_begin;
                        k_local_end = intersection_end - points_begin;

                    }

                    for (long long k = 0; k < k_local_start; k++) {

                        const auto & loc_data = data_not_in_rank_.at(points_begin + k); 
                        
                        std::complex<double> src_dens = density_not_in_rank.at(points_begin + k);
               
                        const long long patch_other = loc_data[6] / patch_size;

                        if (exclude_pts_RP) {

                            auto pre_it = precomputations_data_.find(patch_other);
                            if (pre_it != precomputations_data_.end()) {
                                if (pre_it->second.count(point)) continue;
                            }

                        }

                        std::complex<double> k_val;

                        _kernel(loc_data[0], loc_data[1], loc_data[2], 
                                xp, yp, zp, 
                                loc_data[3], loc_data[4], loc_data[5], 
                                coupling_parameter_, wavenumber_,
                                src_dens, k_val);

                        solution_[i] += k_val;

                    }

                    for (long long k = k_local_start; k < k_local_end; k++) {

                        long long idx = points_begin + k - point_low_;
                        
                        double sx = x_[idx], sy = y_[idx], sz = z_[idx];
                        double nx = normal_x_[idx], ny = normal_y_[idx], nz = normal_z_[idx];
                        std::complex<double> src_dens = density[idx];
                        
                        const long long patch_other = sorting_[idx] / patch_size;

                        if (exclude_pts_RP) {

                            auto pre_it = precomputations_data_.find(patch_other);
                            if (pre_it != precomputations_data_.end()) {
                                if (pre_it->second.count(point)) continue;
                            }

                        }

                        std::complex<double> k_val;

                        _kernel(sx, sy, sz,
                                xp, yp, zp, 
                                nx, ny, nz,
                                coupling_parameter_, wavenumber_,
                                src_dens, k_val);

                        solution_[i] += k_val;                        

                    }

                    for (long long k = k_local_end; k < npoints; k++) {

                        const auto & loc_data = data_not_in_rank_.at(points_begin + k); 
                        
                        std::complex<double> src_dens = density_not_in_rank.at(points_begin + k);
               
                        const long long patch_other = loc_data[6] / patch_size;

                        if (exclude_pts_RP) {

                            auto pre_it = precomputations_data_.find(patch_other);
                            if (pre_it != precomputations_data_.end()) {
                                if (pre_it->second.count(point)) continue;
                            }

                        }

                        std::complex<double> k_val;

                        _kernel(loc_data[0], loc_data[1], loc_data[2], 
                                xp, yp, zp, 
                                loc_data[3], loc_data[4], loc_data[5], 
                                coupling_parameter_, wavenumber_,
                                src_dens, k_val);

                        solution_[i] += k_val;

                    }

                } 

            }

            }

        }

        template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
                               const std::complex<double>, const std::complex<double>, const std::complex<double>, std::complex<double>&),
                  void _factorization(const double, const std::complex<double>, std::complex<double>&)>
        void LevelDEvaluations(const std::vector<std::complex<double>>& density,
                               const std::unordered_map<long long, std::complex<double>>& density_not_in_rank,
                               std::vector<std::complex<double>>& conesegments_current,
                               std::vector<std::complex<double>>& conesegments_prev)
        {
    
            const long long ncones = levels_.back()->GetTotalNumberOfCones();
            const auto & rel_conesegments = levels_.back()->relconesegmentids_;

            #pragma omp parallel
            {

            long long old_cocentered_morton_box = -1;
            long long npoints = 0;
            long long points_begin = 0;
            long long points_end = 0;
            long long k_local_start, k_local_end;

            std::vector<double> t_x(P_), t_y(P_), t_z(P_), t_radii(P_);
            std::vector<std::complex<double>> t_results(P_);

            #pragma omp for schedule(dynamic)
            for (long long i = 0; i < rel_conesegments.size(); i++) {

                const long long rel_cs = rel_conesegments[i];
                const long long cocentered_morton_box = rel_cs / ncones;
                const long long nonrel_conesegment = rel_cs % ncones;

                for (int j = 0; j < P_; j++) {

                    double x, y, z;

                    IPSCHEME_.GetInterpolationPoint(j, x, y, z);

                    t_radii[j] = levels_.back()->Cheb2Radius(nonrel_conesegment, x);

                    levels_.back()->Cheb2Cart(cocentered_morton_box, nonrel_conesegment, x, y, z);
                    
                    t_x[j] = x;
                    t_y[j] = y;
                    t_z[j] = z;

                    t_results[j] = {0.0, 0.0}; 
                    
                }
                
                if (old_cocentered_morton_box != cocentered_morton_box) {

                    const auto & points_data = levels_.back()->mortonbox2discretizationpoints_2_.at(cocentered_morton_box);
                    points_begin = points_data[0];
                    npoints = points_data[1];
                    points_end = points_begin + npoints;

                    const long long intersection_start = std::max(points_begin, point_low_);
                    const long long intersection_end = std::min(points_end, point_up_);

                    if (intersection_start < intersection_end) {
                        
                        k_local_start = intersection_start - points_begin;
                        k_local_end = intersection_end - points_begin;

                    } else {
                        
                        k_local_start = 0;
                        k_local_end = 0;

                    }

                    old_cocentered_morton_box = cocentered_morton_box;

                }       
 
                const long long coefficients_begin = i * P_;

                for (long long k = 0; k < k_local_start; k++) {

                    const auto & loc_data = data_not_in_rank_[points_begin + k]; 
                    
                    std::complex<double> src_dens = density_not_in_rank.at(points_begin + k);
               
                    #pragma omp simd
                    for (int j = 0; j < P_; j++) {

                        std::complex<double> k_val;

                        _kernel(loc_data[0], loc_data[1], loc_data[2], 
                                t_x[j], t_y[j], t_z[j], 
                                loc_data[3], loc_data[4], loc_data[5], 
                                coupling_parameter_, wavenumber_,
                                src_dens, k_val);

                        t_results[j] += k_val;

                    }

                }

                for (long long k = k_local_start; k < k_local_end; k++) {

                    long long idx = points_begin + k - point_low_;
                    
                    double sx = x_[idx], sy = y_[idx], sz = z_[idx];
                    double nx = normal_x_[idx], ny = normal_y_[idx], nz = normal_z_[idx];
                    std::complex<double> src_dens = density[idx];
               
                    #pragma omp simd
                    for (int j = 0; j < P_; j++) {

                        std::complex<double> k_val;

                        _kernel(sx, sy, sz, 
                                t_x[j], t_y[j], t_z[j], 
                                nx, ny, nz, 
                                coupling_parameter_, wavenumber_, 
                                 src_dens, k_val);

                        t_results[j] += k_val;

                    }

                }

                for (long long k = k_local_end; k < npoints; k++) {

                    const auto & loc_data = data_not_in_rank_[points_begin + k]; 
                    
                    std::complex<double> src_dens = density_not_in_rank.at(points_begin + k);
               
                    #pragma omp simd
                    for (int j = 0; j < P_; j++) {

                        std::complex<double> k_val;

                        _kernel(loc_data[0], loc_data[1], loc_data[2], 
                                t_x[j], t_y[j], t_z[j], 
                                loc_data[3], loc_data[4], loc_data[5], 
                                coupling_parameter_, wavenumber_,
                                src_dens, k_val);

                        t_results[j] += k_val;

                    }

                }
                
                for (int j = 0; j < P_; j++) {

                    std::complex<double> fac;
                    _factorization(t_radii[j], wavenumber_, fac);
                    
                    std::complex<double> inv_fac = 1.0 / fac;
                    conesegments_prev[coefficients_begin + j] = t_results[j] * inv_fac;

                }

                IPSCHEME_.GenerateInterpolant(&conesegments_prev[coefficients_begin]);

            }

            }

            std::memcpy(conesegments_current.data(), conesegments_prev.data(), conesegments_current.size() * sizeof(std::complex<double>));

        }

        void CommunicateInterpolation(const int level, std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
                                      std::vector<std::complex<double>>& conesegments_prev) const 
        {

            long long N = conesegments_prev.size();

            std::complex<double>* buffer = conesegments_prev.data();

            MPI_Win win;
            MPI_Win_create(buffer, N * sizeof(std::complex<double>), sizeof(std::complex<double>), MPI_INFO_NULL, mpi_comm_, &win);

            MPI_Win_lock_all(MPI_MODE_NOCHECK, win);

            coeffs.reserve(levels_[level]->flat_send_buffer_interpolation_.size());

            std::vector<MPI_Aint> displacements;
            displacements.reserve(N / comm_size_);

            std::vector<std::complex<double>> recv_buffer;
            recv_buffer.reserve((N / comm_size_) * P_); 
            
            for (int rank = 0; rank < comm_size_; rank++) {

                long long start = levels_[level]->sdispls_interpolation_[rank];
                long long size = levels_[level]->send_counts_interpolation_[rank];

                displacements.resize(size);
                const long long* pos_ptr = &levels_[level]->flat_pos_interpolation_[start];

                #pragma omp simd
                for (long long i = 0; i < size; i++) {
                    
                    displacements[i] = (MPI_Aint)(pos_ptr[i] * sizeof(std::complex<double>));

                }

                MPI_Datatype target_type;
                
                MPI_Type_create_hindexed_block(size, P_, displacements.data(), MPI_DOUBLE_COMPLEX, &target_type);
                MPI_Type_commit(&target_type);

                long long total_elements = size * P_;
                
                recv_buffer.resize(total_elements);

                MPI_Get(recv_buffer.data(), total_elements, MPI_DOUBLE_COMPLEX, rank, 0, 1, target_type, win);

                MPI_Win_flush(rank, win);
                MPI_Type_free(&target_type);            

                const long long* cs_ptr = &levels_[level]->flat_send_buffer_interpolation_[start];
                const std::complex<double>* raw_data = recv_buffer.data();
                
                for (long long i = 0; i < size; i++) {

                    long long cs = cs_ptr[i];
                    
                    const std::complex<double>* seg_start = &raw_data[i * P_];
                    const std::complex<double>* seg_end = seg_start + P_;
                    
                    std::vector<std::complex<double>>& entry = coeffs[cs];
                    entry.assign(seg_start, seg_end);
                }

            }

            MPI_Win_unlock_all(win);

            MPI_Win_free(&win);

        }

        template<void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
                              const std::complex<double>, const std::complex<double>, const std::complex<double>, std::complex<double>&),
                 void _factorization(const double, const std::complex<double>, std::complex<double>&)>
        void Interpolation(const int level, const std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
                           const std::vector<std::complex<double>>& density, const std::unordered_map<long long, std::complex<double>>& density_not_in_rank,
                           std::vector<std::complex<double>>& conesegments_prev,
                           const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
                           const bool exclude_pts_RP)
        {

            const long long ncones = levels_[level]->GetTotalNumberOfCones();
            const auto & split_points = levels_[level]->split_points_relconesegments_;
            const auto & rel_conesegments = levels_[level]->relconesegmentids_;

            const long long N = x.size();

            const long long patch_size = N_PTS_PER_PATCH[0] * N_PTS_PER_PATCH[1];

            #pragma omp parallel
            {

            long long old_morton_box = -1;
            std::vector<long long> morton_cousinboxes;

            #pragma omp for schedule(dynamic)
            for (long long i = 0; i < N; i++) {

                long long point;
                if (exclude_pts_RP) point = sorting_[i];

                const double xp = x[i];
                const double yp = y[i];
                const double zp = z[i];

                const long long mortonbox = levels_[level]->Point2Morton(xp, yp, zp); 

                if (mortonbox != old_morton_box) {

                    morton_cousinboxes = levels_[level]->GetCousins(mortonbox);

                    old_morton_box = mortonbox;

                }

                for (const long long morton_cousinbox : morton_cousinboxes) {

                    const std::array<long long, 2> points_data = levels_[level]->mortonbox2discretizationpoints_2_.at(morton_cousinbox);
                    
                    const long long points_begin = points_data[0];
                    const long long npoints = points_data[1];                          

                    if (USE_ADAPTIVITY && npoints <= MAX_ELEMS_LEAF) {

                        const long long points_end = points_begin + npoints;

                        long long k_local_start = 0, k_local_end = 0;

                        const long long intersection_start = std::max(points_begin, point_low_);
                        const long long intersection_end = std::min(points_end, point_up_);

                        if (intersection_start < intersection_end) {
                            
                            k_local_start = intersection_start - points_begin;
                            k_local_end = intersection_end - points_begin;

                        }

                        for (long long k = 0; k < k_local_start; k++) {

                            const auto & loc_data = data_not_in_rank_[points_begin + k]; 

                            std::complex<double> src_dens = density_not_in_rank.at(points_begin + k);

                            const long long patch_other = loc_data[6] / patch_size;

                            if (exclude_pts_RP) {

                                auto pre_it = precomputations_data_.find(patch_other);
                                if (pre_it != precomputations_data_.end()) {
                                    if (pre_it->second.count(point)) continue;
                                }
    
                            }
    
                            std::complex<double> k_val;
    
                            _kernel(loc_data[0], loc_data[1], loc_data[2], 
                                    xp, yp, zp, 
                                    loc_data[3], loc_data[4], loc_data[5], 
                                    coupling_parameter_, wavenumber_,
                                    src_dens, k_val);
    
                            solution_[i] += k_val;

                        }

                        for (long long k = k_local_start; k < k_local_end; k++) {

                            long long idx = points_begin + k - point_low_;
                        
                            double sx = x_[idx], sy = y_[idx], sz = z_[idx];
                            double nx = normal_x_[idx], ny = normal_y_[idx], nz = normal_z_[idx];
                            std::complex<double> src_dens = density[idx];
                            
                            const long long patch_other = sorting_[idx] / patch_size;

                            if (exclude_pts_RP) {

                                auto pre_it = precomputations_data_.find(patch_other);
                                if (pre_it != precomputations_data_.end()) {
                                    if (pre_it->second.count(point)) continue;
                                }

                            }

                            std::complex<double> k_val;

                            _kernel(sx, sy, sz,
                                    xp, yp, zp, 
                                    nx, ny, nz,
                                    coupling_parameter_, wavenumber_,
                                     src_dens, k_val);

                            solution_[i] += k_val;   

                        }

                        for (long long k = k_local_end; k < npoints; k++) {

                            const auto & loc_data = data_not_in_rank_.at(points_begin + k); 
                        
                            std::complex<double> src_dens = density_not_in_rank.at(points_begin + k);
                
                            const long long patch_other = loc_data[6] / patch_size;

                            if (exclude_pts_RP) {

                                auto pre_it = precomputations_data_.find(patch_other);
                                if (pre_it != precomputations_data_.end()) {
                                    if (pre_it->second.count(point)) continue;
                                }

                            }

                            std::complex<double> k_val;

                            _kernel(loc_data[0], loc_data[1], loc_data[2], 
                                    xp, yp, zp, 
                                    loc_data[3], loc_data[4], loc_data[5], 
                                    coupling_parameter_, wavenumber_,
                                    src_dens, k_val);

                            solution_[i] += k_val;

                        }

                    } else {

                        double locx = xp;
                        double locy = yp;
                        double locz = zp;

                        long long nonrelconesegment = levels_[level]->Cart2Nonrelcone(morton_cousinbox, locx, locy, locz);
                        levels_[level]->Cart2Cheb(morton_cousinbox, locx, locy, locz, nonrelconesegment);
                        const long long relconesegment = morton_cousinbox * ncones + nonrelconesegment;
                        const double radius = levels_[level]->Cheb2Radius(nonrelconesegment, locx);
                        std::complex<double> fac;
                        _factorization(radius, wavenumber_, fac); 

                        const std::complex<double>* vals_ptr = nullptr;
                    
                        if ((relconesegment >= split_points[comm_rank_]) && (relconesegment < split_points[comm_rank_+1])) {
                        
                            auto it = std::lower_bound(rel_conesegments.begin(), rel_conesegments.end(), relconesegment);

                            const long long coeffs_begin_id = std::distance(rel_conesegments.begin(), it) * P_;
                            
                            vals_ptr = &conesegments_prev[coeffs_begin_id];

                        } else {

                            vals_ptr = coeffs.at(relconesegment).data();
                            
                        }   

                        const std::complex<double> tmp = IPSCHEME_.Interpolate(locx, locy, locz, vals_ptr);
                        
                        const double value_real = tmp.real() * fac.real() - tmp.imag() * fac.imag();
                        const double value_imag = tmp.real() * fac.imag() + tmp.imag() * fac.real();

                        solution_[i] += std::complex<double>{value_real, value_imag};

                    }

                }

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
                                    std::vector<std::complex<double>>& conesegments_current) const {

            long long N = conesegments_current.size();

            std::complex<double>* buffer = conesegments_current.data();

            MPI_Win win;
            MPI_Win_create(buffer, N * sizeof(std::complex<double>), sizeof(std::complex<double>), MPI_INFO_NULL, mpi_comm_, &win);

            MPI_Win_lock_all(MPI_MODE_NOCHECK, win);

            coeffs.reserve(levels_[level]->flat_send_buffer_propagation_.size());

            std::vector<MPI_Aint> displacements;
            displacements.reserve(N / comm_size_);

            std::vector<std::complex<double>> recv_buffer;
            recv_buffer.reserve((N / comm_size_) * P_); 
            
            for (int rank = 0; rank < comm_size_; rank++) {

                long long start = levels_[level]->sdispls_propagation_[rank];
                long long size = levels_[level]->send_counts_propagation_[rank];

                displacements.resize(size);
                const long long* pos_ptr = &levels_[level]->flat_pos_propagation_[start];

                #pragma omp simd
                for (long long i = 0; i < size; i++) {
                    
                    displacements[i] = (MPI_Aint)(pos_ptr[i] * sizeof(std::complex<double>));

                }

                MPI_Datatype target_type;
                
                MPI_Type_create_hindexed_block(size, P_, displacements.data(), MPI_DOUBLE_COMPLEX, &target_type);
                MPI_Type_commit(&target_type);

                long long total_elements = size * P_;
                
                recv_buffer.resize(total_elements);

                MPI_Get(recv_buffer.data(), total_elements, MPI_DOUBLE_COMPLEX, rank, 0, 1, target_type, win);

                MPI_Win_flush(rank, win);
                MPI_Type_free(&target_type);            

                const long long* cs_ptr = &levels_[level]->flat_send_buffer_propagation_[start];
                const std::complex<double>* raw_data = recv_buffer.data();
                
                for (long long i = 0; i < size; i++) {

                    long long cs = cs_ptr[i];
                    
                    const std::complex<double>* seg_start = &raw_data[i * P_];
                    const std::complex<double>* seg_end = seg_start + P_;
                    
                    std::vector<std::complex<double>>& entry = coeffs[cs];
                    entry.assign(seg_start, seg_end);
                }

            }

            MPI_Win_unlock_all(win);

            MPI_Win_free(&win);

        }

        template<void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
                              const std::complex<double>, const std::complex<double>, const std::complex<double>, std::complex<double>&),
                 void _factorization(const double, const std::complex<double>, std::complex<double>&)>
        void Propagation(const int level, const std::unordered_map<long long, std::vector<std::complex<double>>>& coeffs,
                            const std::vector<std::complex<double>>& density, const std::unordered_map<long long, std::complex<double>>& density_not_in_rank,
                            std::vector<std::complex<double>>& conesegments_current, const std::vector<std::complex<double>>& conesegments_prev) 
        {

            const long long ncones = levels_[level-1]->GetTotalNumberOfCones();
            const long long ncones_2 = levels_[level]->GetTotalNumberOfCones();
            const auto & rel_conesegments = levels_[level-1]->relconesegmentids_;
            const auto & rel_conesegments_2 = levels_[level]->relconesegmentids_;
            const auto & split_points_2 = levels_[level]->split_points_relconesegments_;

            #pragma omp parallel
            {

            long long old_cocentered_morton_box = -1;
            std::vector<long long> morton_children;    
            std::vector<double> x(P_), y(P_), z(P_), radii(P_);
            std::vector<std::complex<double>> value(P_);                
            
            #pragma omp for schedule(dynamic)
            for (long long i = 0; i < rel_conesegments.size(); i++) {

                const long long rel_cs = rel_conesegments[i];
                const long long cocentered_morton_box = rel_cs / ncones;
                const long long nonrel_conesegment = rel_cs % ncones;

                const long long new_coeffs_begin_id = i * P_; 
                
                if (cocentered_morton_box != old_cocentered_morton_box) {

                    morton_children = GetMortonRelChildren(level-1, cocentered_morton_box);

                    old_cocentered_morton_box = cocentered_morton_box;

                } 

                for (int interpid = 0; interpid < P_; interpid++) {

                    IPSCHEME_.GetInterpolationPoint(interpid, x[interpid], y[interpid], z[interpid]);
                    
                    radii[interpid] = levels_[level-1]->Cheb2Radius(nonrel_conesegment, x[interpid]);

                    levels_[level-1]->Cheb2Cart(cocentered_morton_box, nonrel_conesegment, x[interpid], y[interpid], z[interpid]);

                    value[interpid] = {0.0, 0.0};

                }

                for (long long morton_child : morton_children) {
                    
                    std::array<long long, 2> data_box = levels_[level]->mortonbox2discretizationpoints_2_.at(morton_child);
                    const long long points_begin = data_box[0];
                    const long long npoints = data_box[1];

                    if (USE_ADAPTIVITY && npoints <= MAX_ELEMS_LEAF) {

                        const long long points_end = points_begin + npoints;

                        long long k_local_start = 0, k_local_end = 0;

                        const long long intersection_start = std::max(points_begin, point_low_);
                        const long long intersection_end = std::min(points_end, point_up_);

                        if (intersection_start < intersection_end) {
                            
                            k_local_start = intersection_start - points_begin;
                            k_local_end = intersection_end - points_begin;

                        }

                        for (long long k = 0; k < k_local_start; k++) {

                            const auto & loc_data = data_not_in_rank_[points_begin + k]; 

                            std::complex<double> src_density = density_not_in_rank.at(points_begin + k);
                            
                            #pragma omp simd
                            for (int j = 0; j < P_; j++) {

                                std::complex<double> k_val;

                                _kernel(loc_data[0], loc_data[1], loc_data[2], 
                                        x[j], y[j], z[j], 
                                        loc_data[3], loc_data[4], loc_data[5], 
                                        coupling_parameter_, wavenumber_,
                                        src_density, k_val);

                                value[j] += k_val;

                            }

                        }

                        for (long long k = k_local_start; k < k_local_end; k++) {

                            long long idx = points_begin + k - point_low_;
                            
                            double sx = x_[idx], sy = y_[idx], sz = z_[idx];
                            double nx = normal_x_[idx], ny = normal_y_[idx], nz = normal_z_[idx];
                            std::complex<double> s_dens = density[idx];
                            
                            #pragma omp simd
                            for (int j = 0; j < P_; j++) {

                                std::complex<double> k_val;

                                _kernel(sx, sy, sz, 
                                        x[j], y[j], z[j], 
                                        nx, ny, nz, 
                                        coupling_parameter_, wavenumber_,
                                        s_dens, k_val);

                                value[j] += k_val;

                            }

                        }

                        for (long long k = k_local_end; k < npoints; k++) {

                            const auto & loc_data = data_not_in_rank_[points_begin + k]; 

                            std::complex<double> src_density = density_not_in_rank.at(points_begin + k);
    
                            #pragma omp simd
                            for (int j = 0; j < P_; j++) {

                                std::complex<double> k_val;

                                _kernel(loc_data[0], loc_data[1], loc_data[2], 
                                        x[j], y[j], z[j], 
                                        loc_data[3], loc_data[4], loc_data[5], 
                                        coupling_parameter_, wavenumber_,
                                        src_density, k_val);

                                value[j] += k_val;

                            }

                        }

                    } else {

                        for (int j = 0; j < P_; j++) {

                            double locx = x[j];
                            double locy = y[j];
                            double locz = z[j];

                            long long locnonrel_conesegment = levels_[level]->Cart2Nonrelcone(morton_child, locx, locy, locz);                                
                            levels_[level]->Cart2Cheb(morton_child, locnonrel_conesegment, locx, locy, locz);
                            const long long locrelconesegment = morton_child * ncones_2 + locnonrel_conesegment;
                            const double locradius = levels_[level]->Cheb2Radius(locnonrel_conesegment, locx);
                            std::complex<double> locfac;
                            _factorization(locradius, wavenumber_, locfac);

                            const std::complex<double>* vals_ptr = nullptr;
                    
                            if ((locrelconesegment >= split_points_2[comm_rank_]) && (locrelconesegment < split_points_2[comm_rank_+1])) {
                            
                                auto it = std::lower_bound(rel_conesegments_2.begin(), rel_conesegments_2.end(), locrelconesegment);

                                const long long coeffs_begin_id = std::distance(rel_conesegments_2.begin(), it) * P_;
                                
                                vals_ptr = &conesegments_prev[coeffs_begin_id];

                            } else {

                                vals_ptr = coeffs.at(locrelconesegment).data();
                                
                            }  

                            const std::complex<double> tmp = IPSCHEME_.Interpolate(locx, locy, locz, vals_ptr);

                            const double value_real = tmp.real() * locfac.real() - tmp.imag() * locfac.imag();
                            const double value_imag = tmp.imag() * locfac.real() + tmp.real() * locfac.imag();
                            
                            value[j] += std::complex<double>(value_real, value_imag);

                        }

                    }

                }

                for (int interpid = 0; interpid < P_; interpid++) {

                    std::complex<double> fac;
                    _factorization(radii[interpid], wavenumber_, fac);
                    double norm_sq = std::norm(fac);
                    std::complex<double> result_div_fac = (value[interpid] * std::conj(fac)) / norm_sq;

                    conesegments_current[new_coeffs_begin_id + interpid] = result_div_fac;

                }

                IPSCHEME_.GenerateInterpolant(&conesegments_current[new_coeffs_begin_id]);

            }

            }

        }

    public:

        BoxTree()
        {

        }

        ~BoxTree() 
        {

            for (int iter = 0; iter < levels_.size(); iter++) {

                delete levels_[iter];

            }

        }

        void InitializeObject(typename std::vector<double>::const_iterator x_begin, typename std::vector<double>::const_iterator x_end,
                              typename std::vector<double>::const_iterator y_begin, typename std::vector<double>::const_iterator y_end,
                              typename std::vector<double>::const_iterator z_begin, typename std::vector<double>::const_iterator z_end,
                              typename std::vector<double>::const_iterator nx_begin, typename std::vector<double>::const_iterator nx_end,
                              typename std::vector<double>::const_iterator ny_begin, typename std::vector<double>::const_iterator ny_end,
                              typename std::vector<double>::const_iterator nz_begin, typename std::vector<double>::const_iterator nz_end,
                              const std::vector<long long>& split_points,
                              const std::vector<MPI_Count>& recv_counts,
                              const std::vector<MPI_Aint>& displs,
                              std::complex<double> coupling_parameter, std::complex<double> WAVE_NUMBER, double bbsizeoffset,
                              int Nu_int, int Nv_int, 
                              bool use_adaptivity, bool use_acc, long long max_elems_per_leaf, int n_levels_ifgf, 
                              std::vector<long long>& sorting, MPI_Comm mpi_comm)
        {
            
            x_ = std::vector<double>(x_begin, x_end);
            y_ = std::vector<double>(y_begin, y_end);
            z_ = std::vector<double>(z_begin, z_end);
            
            normal_x_ = std::vector<double>(nx_begin, nx_end);
            normal_y_ = std::vector<double>(ny_begin, ny_end);
            normal_z_ = std::vector<double>(nz_begin, nz_end);

            nlevels_ = n_levels_ifgf;
            wavenumber_ = WAVE_NUMBER;

            mpi_comm_ = mpi_comm;
            N_PTS_PER_PATCH[0] = Nu_int;
            N_PTS_PER_PATCH[1] = Nv_int;
            USE_ADAPTIVITY = use_adaptivity;
            MAX_ELEMS_LEAF = max_elems_per_leaf;
            USE_ACCELERATOR = use_acc;
            BBSIZEOFFSET = bbsizeoffset;


            split_points_orig_ = split_points;
            recv_counts_orig_  = recv_counts;
            displs_orig_       = displs;

            N_loc_orig_ = x_.size();
            coupling_parameter_ = coupling_parameter;

            P_ = (PS * PT * PT);

            Initialize();

            sorting = sorting_;

        }
        
        long long get_box(const double x, const double y, const double z)
        {

            return levels_.back()->Point2Morton(x, y, z);

        }

        void set_precomputations_data(std::unordered_map<long long, std::unordered_set<long long>> data) 
        {
        
            precomputations_data_ = std::move(data);
    
        }
        
        std::vector<long long> get_neighbours_box(const long long i) 
        {

            return levels_.back()->GetNeighboursAll(i);

        }

        std::unordered_map<long long, std::array<long long, 2>> get_mortonbox2discretizationpoints_all()
        {

            std::vector<long long> morton_box_loc;
            std::vector<long long> position_loc;
            std::vector<long long> size_loc;

            morton_box_loc.reserve(levels_.back()->mortonbox2discretizationpoints_.size());
            position_loc.reserve(levels_.back()->mortonbox2discretizationpoints_.size());
            size_loc.reserve(levels_.back()->mortonbox2discretizationpoints_.size());

            for (const auto & [morton_box, values] : levels_.back()->mortonbox2discretizationpoints_) {

                morton_box_loc.push_back(morton_box);
                position_loc.push_back(values[0]);
                size_loc.push_back(values[1]);

            }

            MPI_Count total_size_loc = morton_box_loc.size();

            std::vector<MPI_Count> recv_counts(comm_size_);
            std::vector<MPI_Aint> displs(comm_size_);

            MPI_Allgather(&total_size_loc, 1, MPI_COUNT, &recv_counts[0], 1, MPI_COUNT, mpi_comm_);

            MPI_Count total_size = 0;

            for (int i = 0; i < comm_size_; i++) {

                displs[i] = total_size;
                total_size += recv_counts[i];

            }

            std::vector<long long> morton_box_all(static_cast<size_t>(total_size));
            std::vector<long long> position_all(static_cast<size_t>(total_size));
            std::vector<long long> size_all(static_cast<size_t>(total_size));

            MPI_Allgatherv_c(&morton_box_loc[0], total_size_loc, MPI_LONG_LONG, &morton_box_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
            MPI_Allgatherv_c(&position_loc[0], total_size_loc, MPI_LONG_LONG, &position_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
            MPI_Allgatherv_c(&size_loc[0], total_size_loc, MPI_LONG_LONG, &size_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);

            std::vector<long long>().swap(morton_box_loc);
            std::vector<long long>().swap(position_loc);
            std::vector<long long>().swap(size_loc);

            std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_all;

            for (long long i = 0; i < total_size; i++) {

                const long long morton_box = morton_box_all[i];
                const long long position = position_all[i];
                const long long size = size_all[i];

                auto [it, inserted] = mortonbox2discretizationpoints_all.try_emplace(morton_box, std::array<long long, 2>{position, size});

                if (!inserted) {

                    auto & entry = it->second;

                    entry[0] = std::min(entry[0], position);
                    entry[1] += size;

                }

            } 

            std::vector<long long>().swap(morton_box_all);
            std::vector<long long>().swap(position_all);
            std::vector<long long>().swap(size_all); 

            return mortonbox2discretizationpoints_all;

        }

        template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        const std::complex<double>, const std::complex<double>, const std::complex<double>, std::complex<double>&),
        void _factorization(const double, const std::complex<double>, std::complex<double>&)>
        void Solve(std::vector<std::complex<double>>& density,
                   const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& z,
                   const bool exclude_pts_RP) {

            SolveOrig<_kernel, _factorization>(density, x, y, z, exclude_pts_RP);

        }

        template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        const std::complex<double>, const std::complex<double>, const std::complex<double>, std::complex<double>&),
        void _factorization(const double, const std::complex<double>, std::complex<double>&)>
        void Solve(std::vector<std::complex<double>>& density) {

            Solve<_kernel, _factorization>(density, x_, y_, z_, true);

        }

        template <void _kernel(const double, const double, const double, const double, const double, const double, const double, const double, const double,
        const std::complex<double>, const std::complex<double>, const std::complex<double>, std::complex<double>&),
        void _factorization(const double, const std::complex<double>, std::complex<double>&)>
        void SolveOrig(std::vector<std::complex<double>>& density,
                       const std::vector<double>& x,
                       const std::vector<double>& y,
                       const std::vector<double>& z,
                       const bool exclude_pts_RP) 
        {

            if (nlevels_ < 3) 
                throw std::invalid_argument("IFGF accelerator requires at least 3 levels.");

            const int D = nlevels_ - 1;   
            const long long N = x.size();

            std::unordered_map<long long, std::complex<double>> density_not_in_rank;            
            
            GetDensityNotInRank(density, density_not_in_rank); 
            
            ZeroSolution(N);   

            SingularInteractions<_kernel>(density, density_not_in_rank, x, y, z, exclude_pts_RP);

            long long maxncoeffs = 0;

            for (int level = 2; level < nlevels_; level++) {

                maxncoeffs = std::max<long long>(maxncoeffs, levels_[level]->relconesegmentids_.size());

            }

            std::vector<std::complex<double>> conesegments_current = std::vector<std::complex<double>>(maxncoeffs * P_, {0.0, 0.0});
            std::vector<std::complex<double>> conesegments_prev = std::vector<std::complex<double>>(maxncoeffs * P_, {0.0, 0.0});   
        
            std::unordered_map<long long, std::vector<std::complex<double>>> tmpinterpolationcoeffs;         
            std::unordered_map<long long, std::vector<std::complex<double>>> tmppropagationcoeffs;

            LevelDEvaluations<_kernel, _factorization>(density, density_not_in_rank, conesegments_current, conesegments_prev);

            CommunicatePropagation(D, tmppropagationcoeffs, conesegments_current);

            for (int level = D; level >= 2; level--) {

                CommunicateInterpolation(level, tmpinterpolationcoeffs, conesegments_prev);

                if (level > 2) {
                
                    Propagation<_kernel, _factorization>(level, tmppropagationcoeffs, density, density_not_in_rank, conesegments_current, conesegments_prev);

                    if (level > 3) {

                        CommunicatePropagation(level-1, tmppropagationcoeffs, conesegments_current);          

                    }

                }

                Interpolation<_kernel, _factorization>(level, tmpinterpolationcoeffs, density, density_not_in_rank, conesegments_prev, x, y, z, exclude_pts_RP);

                SwapCoefficients(conesegments_current, conesegments_prev);

            }

            std::vector<std::complex<double>>().swap(conesegments_current);
            std::vector<std::complex<double>>().swap(conesegments_prev);

            tmppropagationcoeffs.clear();
            tmpinterpolationcoeffs.clear();

            density_not_in_rank.clear();

            density = std::move(solution_);          
     
        }

};


#endif