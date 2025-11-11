#include "../solver2.h"

using namespace Eigen;

void Solver::setup_IFGF_choose_levels() {

    int  nlevels_old = nlevels_;
    bool loc_cond_satisfied;
    bool glob_cond_satisfied;

    while (true) {
        create_IFGF_object();
        compute_new_order_points_RP();
        loc_cond_satisfied = check_patch_in_neighbours();

        MPI_Allreduce(&loc_cond_satisfied, &glob_cond_satisfied, 1, MPI_C_BOOL, MPI_LAND, mpi_comm_);

        // Case if the number of levels is too low, then just set it to 3
        if (nlevels_ < 3) {
            nlevels_ = 3;
            break;
         // Case where we go from to many levels to just the right amount
        } else if(glob_cond_satisfied && nlevels_old > nlevels_) {

             break;

        // Case where on the first try it worked, so we try to go deeper
        } else if (glob_cond_satisfied && nlevels_old <= nlevels_) {
            nlevels_old = nlevels_;
            nlevels_ += 1;

        // Case where the first try it didn't work, so we move down a level
        } else if (!glob_cond_satisfied) {
            nlevels_old = nlevels_;
            nlevels_ -= 1;
        } 

    }


}

void Solver::create_IFGF_object() 
{

    std::vector<double> x, y, z;
    std::vector<double> normal_x, normal_y, normal_z;

    x.insert(x.end(), disc_points_x_all_.begin() + point_low_, disc_points_x_all_.begin() + point_up_);
    y.insert(y.end(), disc_points_y_all_.begin() + point_low_, disc_points_y_all_.begin() + point_up_);
    z.insert(z.end(), disc_points_z_all_.begin() + point_low_, disc_points_z_all_.begin() + point_up_);

    normal_x.insert(normal_x.end(), norm_points_x_all_.begin() + point_low_, norm_points_x_all_.begin() + point_up_);
    normal_y.insert(normal_y.end(), norm_points_y_all_.begin() + point_low_, norm_points_y_all_.begin() + point_up_);
    normal_z.insert(normal_z.end(), norm_points_z_all_.begin() + point_low_, norm_points_z_all_.begin() + point_up_);

    boxes_.InitializeObject(x, y, z, normal_x, normal_y, normal_z, 
                            split_points_2_, recv_counts_2_, displs_2_,
                            coupling_parameter_,
                            new_order_points_IFGF_,
                            nlevels_, wavenumber_,
                            mpi_comm_);

    std::vector<long long> point_precomputations_all;
    std::vector<long long> patch_num_precomputations_all;

    int size_loc = point_precomputations_.size();

    std::vector<int> recv_counts(comm_size_);
    std::vector<int> displs(comm_size_, 0);

    MPI_Allgather(&size_loc, 1, MPI_INT, &recv_counts[0], 1, MPI_INT, mpi_comm_);

    int size_all = recv_counts[0];
    for (int i = 1; i < comm_size_; i++) {
        size_all += recv_counts[i];
        displs[i] = displs[i-1] + recv_counts[i-1];
    }

    point_precomputations_all.resize(size_all);
    patch_num_precomputations_all.resize(size_all);

    MPI_Allgatherv(&point_precomputations_[0], size_loc, MPI_LONG_LONG, &point_precomputations_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);
    MPI_Allgatherv(&patch_num_precomputations_[0], size_loc, MPI_LONG_LONG, &patch_num_precomputations_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);

    std::unordered_map<long long, std::unordered_set<long long>> precomputations_data_all;

    for (long long i = 0; i < size_all; i++) {

        const long long point = point_precomputations_all[i];
        const long long patch = patch_num_precomputations_all[i];

        precomputations_data_all[point].insert(patch);

    }

    std::vector<long long>().swap(point_precomputations_all);
    std::vector<long long>().swap(patch_num_precomputations_all);

    std::unordered_map<long long, std::unordered_set<long long>> precomputations_data_loc;

    for (long long i = point_low_; i < point_up_; i++) {

        const long long point = new_order_points_IFGF_[i];

        precomputations_data_loc[point] = precomputations_data_all[point];

    }

    precomputations_data_all.clear();

    boxes_.set_precomputations_data(precomputations_data_loc);

    precomputations_data_loc.clear();

    boxes_.set_n_pts_per_patch(Nu_int_, Nv_int_);

}

void Solver::compute_new_order_points_RP() {

    std::vector<std::vector<long long>> points_not_in_rank(comm_size_), order_points_not_in_rank(comm_size_);

    new_order_points_RP_ = std::vector<long long>(point_up_-point_low_);

    for (long long i = point_low_; i < point_up_; i++) {

        const long long point = new_order_points_IFGF_[i];

        int rank;

        for (int j = 0; j < comm_size_+1; j++) {

            if ((point >= split_points_2_[j]) && (point < split_points_2_[j+1])) {

                rank = j;
                break;

            }

        }

        if (comm_rank_ == rank) {

            new_order_points_RP_[point-point_low_] = i;

        } else {

            points_not_in_rank[rank].push_back(point);
            order_points_not_in_rank[rank].push_back(i);

        }

    }

    for (int rank = 0; rank < comm_size_; rank++) {

        std::vector<long long> points = points_not_in_rank[rank];
        std::vector<long long> orders = order_points_not_in_rank[rank];

        std::vector<long long> points_all;
        std::vector<long long> orders_all;

        int size_points = points.size();
        int size_total;
        
        std::vector<int> recv_counts(comm_size_);
        std::vector<int> displs(comm_size_, 0);

        MPI_Gather(&size_points, 1, MPI_INT, &recv_counts[0], 1, MPI_INT, rank, mpi_comm_);

        if (rank == comm_rank_) {

            size_total = recv_counts[0];

            for (int i = 1; i < comm_size_; i++) {

                size_total += recv_counts[i];

                displs[i] = displs[i-1] + recv_counts[i-1];

            }

            points_all.resize(size_total);
            orders_all.resize(size_total);

        }

        MPI_Gatherv(&points[0], size_points, MPI_LONG_LONG, &points_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, rank, mpi_comm_);
        MPI_Gatherv(&orders[0], size_points, MPI_LONG_LONG, &orders_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, rank, mpi_comm_);

        std::vector<long long>().swap(points);
        std::vector<long long>().swap(orders);

        if (rank == comm_rank_) {

            for (long long i = 0; i < size_total; i++) {

                new_order_points_RP_[points_all[i]-point_low_] = orders_all[i];

            }

        }

        std::vector<long long>().swap(orders_all);
        std::vector<long long>().swap(points_all);

    }

}

bool Solver::check_patch_in_neighbours() {

    // Compute point to box

    std::vector<long long> point_to_box(point_up_-point_low_);

    for (long long i = 0; i < point_up_-point_low_; i++) {

        point_to_box[i] = boxes_.get_box(new_order_points_RP_[i]);

    }

    // Compute patch to box

    std::unordered_map<long long, std::set<long long>> patch_to_box;

    std::vector<long long> patches_loc;
    std::vector<long long> boxes_loc;
    std::vector<long long> total_boxes_loc;

    for (long long patch_num = 0; patch_num < patch_up_-patch_low_; patch_num++) {

        std::unordered_set<long long> boxes(point_to_box.begin() + patch_num * Nu_int_*Nv_int_, point_to_box.begin() + (patch_num+1) * Nu_int_*Nv_int_);

        patches_loc.push_back(patch_num + patch_low_);
        total_boxes_loc.push_back(boxes.size());
        boxes_loc.insert(boxes_loc.end(), boxes.begin(), boxes.end());

    }

    int size_1_loc = patches_loc.size();
    int size_2_loc = boxes_loc.size();

    std::vector<int> recv_counts_1(comm_size_);
    std::vector<int> recv_counts_2(comm_size_);
    std::vector<int> displs_1(comm_size_, 0);
    std::vector<int> displs_2(comm_size_, 0);

    MPI_Allgather(&size_1_loc, 1, MPI_INT, &recv_counts_1[0], 1, MPI_INT, mpi_comm_);
    MPI_Allgather(&size_2_loc, 1, MPI_INT, &recv_counts_2[0], 1, MPI_INT, mpi_comm_);

    int size_1_all = recv_counts_1[0];
    int size_2_all = recv_counts_2[0];

    for (int i = 1; i < comm_size_; i++) {

        size_1_all += recv_counts_1[i];
        size_2_all += recv_counts_2[i];

        displs_1[i] = displs_1[i-1] + recv_counts_1[i-1];
        displs_2[i] = displs_2[i-1] + recv_counts_2[i-1];

    }

    std::vector<long long> patches_all(size_1_all);
    std::vector<long long> total_boxes_all(size_1_all);
    std::vector<long long> boxes_all(size_2_all);

    MPI_Allgatherv(&patches_loc[0], size_1_loc, MPI_LONG_LONG, &patches_all[0], &recv_counts_1[0], &displs_1[0], MPI_LONG_LONG, mpi_comm_);
    MPI_Allgatherv(&total_boxes_loc[0], size_1_loc, MPI_LONG_LONG, &total_boxes_all[0], &recv_counts_1[0], &displs_1[0], MPI_LONG_LONG, mpi_comm_);
    MPI_Allgatherv(&boxes_loc[0], size_2_loc, MPI_LONG_LONG, &boxes_all[0], &recv_counts_2[0], &displs_2[0], MPI_LONG_LONG, mpi_comm_);

    std::vector<long long>().swap(patches_loc);
    std::vector<long long>().swap(total_boxes_loc);
    std::vector<long long>().swap(boxes_loc);

    long long counter = 0;

    for (long long i = 0; i < size_1_all; i++) {

        long long patch = patches_all[i];
        long long size_boxes = total_boxes_all[i];

        patch_to_box[patch].insert(boxes_all.begin() + counter, boxes_all.begin() + counter + size_boxes);

        counter += size_boxes;

    }

    std::vector<long long>().swap(patches_all);
    std::vector<long long>().swap(total_boxes_all);
    std::vector<long long>().swap(boxes_all);
    
    // Compute relevant morton boxes

    std::unordered_set<long long> mortonidofrelboxes;
    boxes_.get_mortonidofrelboxes(mortonidofrelboxes);

    std::vector<long long> mortonidofrelboxes_loc(mortonidofrelboxes.begin(), mortonidofrelboxes.end());
    mortonidofrelboxes.clear();

    int size_loc = mortonidofrelboxes_loc.size();

    std::vector<int> recv_counts(comm_size_);
    std::vector<int> displs(comm_size_, 0);

    MPI_Allgather(&size_loc, 1, MPI_INT, &recv_counts[0], 1, MPI_INT, mpi_comm_);

    int size_all = recv_counts[0];

    for (int i = 1; i < comm_size_; i++) {

        displs[i] = displs[i-1] + recv_counts[i-1];

        size_all += recv_counts[i];

    }

    std::vector<long long> mortonidofrelboxes_all(size_all);

    MPI_Allgatherv(&mortonidofrelboxes_loc[0], size_loc, MPI_LONG_LONG, &mortonidofrelboxes_all[0], &recv_counts[0], &displs[0], MPI_LONG_LONG, mpi_comm_);

    std::vector<long long>().swap(mortonidofrelboxes_loc);

    for (long long i = 0; i < size_all; i++) {

        mortonidofrelboxes.insert(mortonidofrelboxes_all[i]);

    }

    std::vector<long long>().swap(mortonidofrelboxes_all);

    // Check if patch is in neighbour union

    long long actual_npoint = -1;
    std::vector<long long> neighbours_box_npoint;

    for (long long i = 0; i < point_precomputations_.size(); i++) {

        const long long npoint = point_precomputations_[i];

        if (actual_npoint != npoint) {
        
            const long long box_npoint = point_to_box[npoint - point_low_];
            neighbours_box_npoint = boxes_.get_neighbours_box(box_npoint);

            std::vector<long long> relevant_neighbours;
            for (long long j = 0; j < neighbours_box_npoint.size(); j++) {

                if (mortonidofrelboxes.count(neighbours_box_npoint[j]) != 0) {
                    relevant_neighbours.push_back(neighbours_box_npoint[j]);
                }

            }

            neighbours_box_npoint = relevant_neighbours;

            actual_npoint = npoint;

        }

        const long long patch_num_singular = patch_num_precomputations_[i];
        const std::set<long long> boxes_patch_num_singular = patch_to_box[patch_num_singular];

        if (!(std::includes(neighbours_box_npoint.begin(), neighbours_box_npoint.end(), boxes_patch_num_singular.begin(), boxes_patch_num_singular.end()))) {

            //std::cout << "Point " << npoint << " and singular / near singular patches are not contained in neighbor union.\n";
            // std::exit(0);
            return false;

        }

    }

    return true;

}



void Solver::initialize_indexes_HO() {

    boxes_.InitializeIndexes(
    [&](double a1,double a2,double a3,double a4,double a5,double a6,
        double a7,double a8,double a9,std::complex<double> a10,
        std::complex<double> c1,std::complex<double> c2,std::complex<double>& out) {
        Solver::fct_4(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,c1,c2,out);
    }
    );

    // boxes_.InitializeIndexes<fct_4>();

}