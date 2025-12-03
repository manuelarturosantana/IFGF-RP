#include "../solver2.h"

void Solver::create_IFGF_object() 
{

    boxes_.InitializeObject(disc_points_x_all_.begin(), disc_points_x_all_.end(),
                            disc_points_y_all_.begin(), disc_points_y_all_.end(),
                            disc_points_z_all_.begin(), disc_points_z_all_.end(),
                            norm_points_x_all_.begin(), norm_points_x_all_.end(),
                            norm_points_y_all_.begin(), norm_points_y_all_.end(),
                            norm_points_z_all_.begin(), norm_points_z_all_.end(),
                            split_points_2_, recv_counts_2_, displs_2_,
                            coupling_parameter_, WAVE_NUMBER, EQUATION_FORMULATION,
                            new_order_points_IFGF_, mpi_comm_);
}

void Solver::set_precomputations_data_IFGF()
{

    std::unordered_map<long long, std::unordered_set<long long>> patch_to_sing_point;
    std::vector<std::vector<long long>> points_not_in_rank(comm_size_);

    for (long long i = 0; i < point_up_ - point_low_; i++) {

        const long long point = new_order_points_IFGF_[i];

        auto it = std::upper_bound(split_points_2_.begin(), split_points_2_.end(), point);

        int rank = static_cast<int>(std::distance(split_points_2_.begin(), it)) - 1;

        if (rank < 0) {
        
            rank = 0;
        
        } else if (rank >= comm_size_) {
        
            rank = comm_size_ - 1; 
        
        }

        if (comm_rank_ == rank) {

            const auto point_begin = std::lower_bound(point_precomputations_.begin(), point_precomputations_.end(), point);
            const auto point_end = std::upper_bound(point_precomputations_.begin(), point_precomputations_.end(), point);

            const size_t start_idx = std::distance(point_precomputations_.begin(), point_begin);
            const size_t end_idx = std::distance(point_precomputations_.begin(), point_end);

            for (size_t j = start_idx; j < end_idx; ++j) {
            
                patch_to_sing_point[patch_num_precomputations_[j]].insert(point);
    
            }

        } else {

            points_not_in_rank[rank].push_back(point);

        }

    }

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

    points_not_in_rank.clear();

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

    std::vector<std::vector<MPI_Count>> patch_counts_to_send_back(comm_size_);
    std::vector<std::vector<long long>> patches_data_to_send_back(comm_size_);

    MPI_Aint current_flat_recv_points_idx = 0;

    for (int sender_rank = 0; sender_rank < comm_size_; ++sender_rank) {

        MPI_Count num_points_from_this_sender = recv_counts_points[sender_rank];

        for (MPI_Count i = 0; i < num_points_from_this_sender; ++i) {
    
            const long long point = flat_recv_points_buffer[current_flat_recv_points_idx + i];

            const auto point_begin = std::lower_bound(point_precomputations_.begin(), point_precomputations_.end(), point);
            const auto point_end = std::upper_bound(point_precomputations_.begin(), point_precomputations_.end(), point);

            const size_t start_idx = std::distance(point_precomputations_.begin(), point_begin);
            const size_t end_idx = std::distance(point_precomputations_.begin(), point_end);

            MPI_Count num_patches_for_this_point = static_cast<MPI_Count>(end_idx - start_idx);

            patch_counts_to_send_back[sender_rank].push_back(num_patches_for_this_point);
    
            for (size_t j = start_idx; j < end_idx; ++j) {
        
                patches_data_to_send_back[sender_rank].push_back(patch_num_precomputations_[j]);
        
            }

        }

        current_flat_recv_points_idx += num_points_from_this_sender;

    }            

    flat_recv_points_buffer.clear();

    std::vector<MPI_Count> send_counts_patch_counts(comm_size_);
    std::vector<MPI_Aint> sdispls_patch_counts(comm_size_);
    MPI_Count total_send_patch_counts_flat = 0;

    for (int i = 0; i < comm_size_; ++i) {

        send_counts_patch_counts[i] = static_cast<MPI_Count>(patch_counts_to_send_back[i].size());
        sdispls_patch_counts[i] = total_send_patch_counts_flat;
        total_send_patch_counts_flat += send_counts_patch_counts[i];

    }

    std::vector<MPI_Count> flat_send_patch_counts_buffer(static_cast<size_t>(total_send_patch_counts_flat));

    MPI_Count current_flat_send_pos = 0;

    for (int i = 0; i < comm_size_; ++i) {

        std::copy(patch_counts_to_send_back[i].begin(), patch_counts_to_send_back[i].end(),
                flat_send_patch_counts_buffer.begin() + current_flat_send_pos);
        current_flat_send_pos += send_counts_patch_counts[i];

    }

    patch_counts_to_send_back.clear();

    std::vector<MPI_Count> recv_counts_patch_counts(comm_size_);

    MPI_Alltoall(send_counts_patch_counts.data(), 1, MPI_COUNT,
                    recv_counts_patch_counts.data(), 1, MPI_COUNT, mpi_comm_);

    MPI_Count total_recv_patch_counts_flat = 0;
    std::vector<MPI_Aint> rdispls_patch_counts(comm_size_);

    for (int i = 0; i < comm_size_; ++i) {

        rdispls_patch_counts[i] = total_recv_patch_counts_flat;
        total_recv_patch_counts_flat += recv_counts_patch_counts[i];

    }

    std::vector<MPI_Count> flat_recv_patch_counts_buffer(static_cast<size_t>(total_recv_patch_counts_flat));

    MPI_Alltoallv_c(flat_send_patch_counts_buffer.data(), send_counts_patch_counts.data(), sdispls_patch_counts.data(), MPI_COUNT,
                    flat_recv_patch_counts_buffer.data(), recv_counts_patch_counts.data(), rdispls_patch_counts.data(), MPI_COUNT,
                    mpi_comm_);

    flat_send_patch_counts_buffer.clear();

    std::vector<MPI_Count> send_counts_patches_data(comm_size_);
    std::vector<MPI_Aint> sdispls_patches_data(comm_size_);
    MPI_Count total_send_patches_data_flat = 0;

    for (int i = 0; i < comm_size_; ++i) {

        send_counts_patches_data[i] = static_cast<MPI_Count>(patches_data_to_send_back[i].size());
        sdispls_patches_data[i] = total_send_patches_data_flat;
        total_send_patches_data_flat += send_counts_patches_data[i];

    }

    std::vector<long long> flat_send_patches_data_buffer(static_cast<size_t>(total_send_patches_data_flat));
    current_flat_send_pos = 0;

    for (int i = 0; i < comm_size_; ++i) {

        std::copy(patches_data_to_send_back[i].begin(), patches_data_to_send_back[i].end(),
                flat_send_patches_data_buffer.begin() + current_flat_send_pos);
        current_flat_send_pos += send_counts_patches_data[i];

    }

    patches_data_to_send_back.clear();

    std::vector<MPI_Count> recv_counts_patches_data(comm_size_);

    MPI_Alltoall(send_counts_patches_data.data(), 1, MPI_COUNT,
                    recv_counts_patches_data.data(), 1, MPI_COUNT, mpi_comm_);

    MPI_Count total_recv_patches_data_flat = 0;
    std::vector<MPI_Aint> rdispls_patches_data(comm_size_);

    for (int i = 0; i < comm_size_; ++i) {

        rdispls_patches_data[i] = total_recv_patches_data_flat;
        total_recv_patches_data_flat += recv_counts_patches_data[i];

    }

    std::vector<long long> flat_recv_patches_data_buffer(static_cast<size_t>(total_recv_patches_data_flat));

    MPI_Alltoallv_c(flat_send_patches_data_buffer.data(), send_counts_patches_data.data(), sdispls_patches_data.data(), MPI_LONG_LONG,
                    flat_recv_patches_data_buffer.data(), recv_counts_patches_data.data(), rdispls_patches_data.data(), MPI_LONG_LONG,
                    mpi_comm_);

    flat_send_patches_data_buffer.clear();

    MPI_Count count = 0;

    for (MPI_Count i = 0; i < total_send_points; i++) {

        const long long point = flat_send_points_buffer[i];
        const long long size = flat_recv_patch_counts_buffer[i];

        for (long long j = 0; j < size; j++) {

            const long long patch = flat_recv_patches_data_buffer[count];

            patch_to_sing_point[patch].insert(point);

            count++;

        }

    }

    flat_send_points_buffer.clear();
    flat_recv_patch_counts_buffer.clear();
    flat_recv_patches_data_buffer.clear();

    boxes_.set_precomputations_data(std::move(patch_to_sing_point));

}

void Solver::compute_new_order_points_RP() 
{

    std::vector<std::vector<long long>> points_not_in_rank(comm_size_), order_points_not_in_rank(comm_size_);

    new_order_points_RP_ = std::vector<long long>(point_up_-point_low_);

    for (long long i = point_low_; i < point_up_; i++) {

        const long long point = new_order_points_IFGF_[i - point_low_];

        auto it = std::upper_bound(split_points_2_.begin(), split_points_2_.end(), point);

        int rank = std::clamp<int>(static_cast<int>(std::distance(split_points_2_.begin(), it)) - 1, 0, comm_size_-1);

        if (comm_rank_ == rank) {

            new_order_points_RP_[point - point_low_] = i;

        } else {

            points_not_in_rank[rank].push_back(point);
            order_points_not_in_rank[rank].push_back(i);

        }

    }

    std::vector<MPI_Count> send_counts_points(comm_size_);
    std::vector<MPI_Aint> sdispls_points(comm_size_);
    MPI_Count total_send_points = 0;    

    for (int i = 0; i < comm_size_; ++i) {

        send_counts_points[i] = static_cast<MPI_Count>(points_not_in_rank[i].size());
        sdispls_points[i] = total_send_points;
        total_send_points += send_counts_points[i];

    }

    std::vector<long long> flat_send_points_buffer(static_cast<size_t>(total_send_points));
    std::vector<long long> flat_send_orders_buffer(static_cast<size_t>(total_send_points));

    MPI_Aint current_point_pos = 0;

    for (int i = 0; i < comm_size_; ++i) {

        std::copy(points_not_in_rank[i].begin(), points_not_in_rank[i].end(),
                    flat_send_points_buffer.begin() + current_point_pos);
        std::copy(order_points_not_in_rank[i].begin(), order_points_not_in_rank[i].end(),
                    flat_send_orders_buffer.begin() + current_point_pos);
        current_point_pos += send_counts_points[i];

    }

    points_not_in_rank.clear();
    order_points_not_in_rank.clear();

    std::vector<MPI_Count> recv_counts_points(comm_size_);

    MPI_Alltoall(send_counts_points.data(), 1, MPI_COUNT, recv_counts_points.data(), 1, MPI_COUNT, mpi_comm_);

    std::vector<MPI_Aint> rdispls_points(comm_size_);
    MPI_Count total_recv_points = 0;

    for (int i = 0; i < comm_size_; ++i) {

        rdispls_points[i] = total_recv_points;
        total_recv_points += recv_counts_points[i];

    }

    std::vector<long long> flat_recv_points_buffer(static_cast<size_t>(total_recv_points));
    std::vector<long long> flat_recv_orders_buffer(static_cast<size_t>(total_recv_points));

    MPI_Alltoallv_c(flat_send_points_buffer.data(), send_counts_points.data(), sdispls_points.data(), MPI_LONG_LONG,
                    flat_recv_points_buffer.data(), recv_counts_points.data(), rdispls_points.data(), MPI_LONG_LONG,
                    mpi_comm_);

    MPI_Alltoallv_c(flat_send_orders_buffer.data(), send_counts_points.data(), sdispls_points.data(), MPI_LONG_LONG, 
                    flat_recv_orders_buffer.data(), recv_counts_points.data(), rdispls_points.data(), MPI_LONG_LONG, 
                    mpi_comm_);

    flat_send_points_buffer.clear();
    flat_send_orders_buffer.clear();

    for (long long i = 0; i < total_recv_points; ++i) {

        const long long point_global_id = flat_recv_points_buffer[i];
        const long long original_global_index = flat_recv_orders_buffer[i];

        new_order_points_RP_[point_global_id - point_low_] = original_global_index;

    }

    flat_recv_points_buffer.clear();
    flat_recv_orders_buffer.clear();

}

void Solver::check_patch_in_neighbours() 
{

    // Compute point to box

    std::vector<long long> point_to_box(point_up_-point_low_);

    #pragma omp parallel for
    for (long long i = 0; i < point_up_-point_low_; i++) {

        point_to_box[i] = boxes_.get_box(disc_points_x_all_[i], disc_points_y_all_[i], disc_points_z_all_[i]);

    }

    // Compute patch to box

    std::vector<std::vector<long long>> patch_to_box(patch_up_-patch_low_);           

    for (long long patch_num = 0; patch_num < patch_up_-patch_low_; patch_num++) {

        auto &vec = patch_to_box[patch_num];

        vec.assign(point_to_box.begin() + patch_num * Nu_int_*Nv_int_, point_to_box.begin() + (patch_num + 1) * Nu_int_*Nv_int_);

        std::sort(vec.begin(), vec.end());
        auto last_unique = std::unique(vec.begin(), vec.end());
        vec.erase(last_unique, vec.end());

    }

    // Send patch_to_box map to all ranks

    std::vector<long long> patches_loc;
    std::vector<long long> size_loc;
    std::vector<long long> boxes_loc;

    patches_loc.reserve(patch_up_-patch_low_);
    size_loc.reserve(patch_up_-patch_low_);
    
    for (long long patch_num = 0; patch_num < patch_up_-patch_low_; patch_num++) {

        patches_loc.push_back(patch_num + patch_low_);
        size_loc.push_back(patch_to_box[patch_num].size());
        boxes_loc.insert(boxes_loc.end(), patch_to_box[patch_num].begin(), patch_to_box[patch_num].end());

    }

    patch_to_box.clear();

    MPI_Count total_1 = patches_loc.size();
    MPI_Count total_2 = boxes_loc.size();

    std::vector<MPI_Count> recv_counts_1(comm_size_);
    std::vector<MPI_Count> recv_counts_2(comm_size_);

    MPI_Allgather(&total_1, 1, MPI_COUNT, &recv_counts_1[0], 1, MPI_COUNT, mpi_comm_);
    MPI_Allgather(&total_2, 1, MPI_COUNT, &recv_counts_2[0], 1, MPI_COUNT, mpi_comm_);

    std::vector<MPI_Aint> displs_1(comm_size_, 0);
    std::vector<MPI_Aint> displs_2(comm_size_, 0);

    MPI_Count total_1_all = 0;
    MPI_Count total_2_all = 0;

    for (int i = 0; i < comm_size_; i++) {

        displs_1[i] = total_1_all;
        displs_2[i] = total_2_all;

        total_1_all += recv_counts_1[i];
        total_2_all += recv_counts_2[i];

    }

    std::vector<long long> patches_all(total_1_all);
    std::vector<long long> size_all(total_1_all);
    std::vector<long long> boxes_all(total_2_all);

    MPI_Allgatherv_c(&patches_loc[0], total_1, MPI_LONG_LONG, &patches_all[0], &recv_counts_1[0], &displs_1[0], MPI_LONG_LONG, mpi_comm_);
    MPI_Allgatherv_c(&size_loc[0], total_1, MPI_LONG_LONG, &size_all[0], &recv_counts_1[0], &displs_1[0], MPI_LONG_LONG, mpi_comm_);
    MPI_Allgatherv_c(&boxes_loc[0], total_2, MPI_LONG_LONG, &boxes_all[0], &recv_counts_2[0], &displs_2[0], MPI_LONG_LONG, mpi_comm_);

    std::vector<long long>().swap(patches_loc);
    std::vector<long long>().swap(size_loc);
    std::vector<long long>().swap(boxes_loc);

    std::vector<std::vector<long long>> patch_to_box_all(Q_ * Qx_*Qy_);   

    long long counter = 0;

    for (long long i = 0; i < Q_ * Qx_*Qy_; i++) {

        long long patch = patches_all[i];
        long long size_boxes = size_all[i];

        auto &vec = patch_to_box_all[i];

        vec.assign(boxes_all.begin() + counter, boxes_all.begin() + counter + size_boxes);

        counter += size_boxes;

    }

    std::vector<long long>().swap(patches_all);
    std::vector<long long>().swap(size_all);
    std::vector<long long>().swap(boxes_all);

    // Compute relevant morton boxes

    std::unordered_set<long long> mortonidofrelboxes;

    for (long long i = 0; i < Q_ * Qx_*Qy_; i++) {

        mortonidofrelboxes.insert(patch_to_box_all[i].begin(), patch_to_box_all[i].end());

    }    

    // Get size boxes level D

    std::unordered_map<long long, std::array<long long, 2>> mortonbox2discretizationpoints_all;

    if (USE_ADAPTIVITY) {

        mortonbox2discretizationpoints_all = boxes_.get_mortonbox2discretizationpoints_all();

    }

    // Check if patch is in neighbour union

    #pragma omp parallel
    {

    long long actual_npoint = -1;
    std::vector<long long> neighbours_box_npoint;

    #pragma omp for
    for (long long i = 0; i < point_precomputations_.size(); i++) {

        const long long npoint = point_precomputations_[i];

        if (actual_npoint != npoint) {
        
            const long long box_npoint = point_to_box[npoint - point_low_];
            neighbours_box_npoint = boxes_.get_neighbours_box(box_npoint);

            std::vector<long long> relevant_neighbours;
            for (const auto & elem : neighbours_box_npoint) {

                if (mortonidofrelboxes.find(elem) != mortonidofrelboxes.end()) {

                    relevant_neighbours.push_back(elem);

                }

            }

            neighbours_box_npoint = std::move(relevant_neighbours);

            actual_npoint = npoint;

        }

        const long long patch_num_singular = patch_num_precomputations_[i];
        const std::vector<long long> boxes_patch_num_singular = patch_to_box_all[patch_num_singular];

        for (const auto & box : boxes_patch_num_singular) {

            bool is_in_adaptive_boxes = true;
            
            if (USE_ADAPTIVITY) {

                is_in_adaptive_boxes = (mortonbox2discretizationpoints_all[box][1] <= MAX_ELEMS_LEAF);

            }                    
            
            bool is_in_neighbour_union = (std::find(neighbours_box_npoint.begin(), neighbours_box_npoint.end(), box) != neighbours_box_npoint.end());

            if (!(is_in_neighbour_union || is_in_adaptive_boxes)) {

                std::cout << "Point " << npoint << " and singular / near singular patches are not contained in neighbor union.\n";
                std::exit(0);

            }

        }

    }

    }

}


