#ifndef BIFROST_HELSITIGS_TCC
#define BIFROST_HELSITIGS_TCC

extern "C" void helsitigs_initialise(size_t threads);

extern "C" void helsitigs_initialise_graph(size_t unitig_amount);

extern "C" void helsitigs_merge_nodes(size_t unitig_a, bool strand_a, size_t unitig_b, bool strand_b);

extern "C" void helsitigs_build_graph(const size_t* unitig_weights);

extern "C" void helsitigs_compute_tigs(size_t tig_algorithm, size_t threads, size_t k, const char* matching_file_prefix, ptrdiff_t* tigs_edge_out, size_t* tigs_insert_out, size_t* tigs_out_limits);



template<typename U, typename G>
bool CompactedDBG<U, G>::convert_tigs(CompactedDBG<U, G>* dbg, const Tigs tigs, const size_t nb_threads, const string& matching_file_prefix) {
    /*cout << "convert_tigs(dbg, tigs, nb_threads = " << nb_threads << ")" << endl;

    cout << "given unitigs:" << endl;
    vector<string> unitigs;
    for (const auto unitig : *dbg) {
        const string seq(unitig.referenceUnitigToString());
        cout << seq << endl;
        unitigs.push_back(seq);
    }*/
    const auto start_send = std::chrono::high_resolution_clock::now();
    cout << "CompactedDBG::convert_tigs(): start" << endl;

    helsitigs_initialise(nb_threads);
    cout << "dbg->size() = " << dbg->size() << endl;
    cout << "dbg->v_unitigs.size() = " << dbg->v_unitigs.size() << endl;
    cout << "dbg->km_unitigs.size() = " << dbg->km_unitigs.size() << endl;
    cout << "dbg->h_kmers_ccov.size() = " << dbg->h_kmers_ccov.size() << endl;
    helsitigs_initialise_graph(dbg->size());
    vector<size_t> unitig_weights;
    auto h_kmer_ccov_ranks = dbg->h_kmers_ccov.compute_rank_array();
    auto h_kmer_ccov_orders = dbg->h_kmers_ccov.compute_order_array();
    /*cout << "h_kmer_ccov_ranks(" << h_kmer_ccov_ranks.size() << ")[";
    for (const auto& r : h_kmer_ccov_ranks) {
        cout << ", " << r;
    }
    cout << "]" << endl;
    cout << "h_kmer_ccov_orders(" << h_kmer_ccov_orders.size() << ")[";
    for (const auto& r : h_kmer_ccov_orders) {
        cout << ", " << r;
    }
    cout << "]" << endl;*/

    const auto start_merge_nodes = std::chrono::high_resolution_clock::now();

    struct MergeCall {size_t u1_index; bool u1_strand; size_t u2_index; bool u2_strand;};

    auto do_merge_calls = [](const vector<MergeCall>& merge_calls) {
        const auto start = std::chrono::high_resolution_clock::now();
        for (const auto& merge_call : merge_calls) {
            helsitigs_merge_nodes(merge_call.u1_index, merge_call.u1_strand, merge_call.u2_index, merge_call.u2_strand);
        }
        const auto stop = std::chrono::high_resolution_clock::now();
        return ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()) / 1e6;
    };

    vector<thread> topology_workers;
    mutex rust_helsitigs_mutex;
    size_t total_merge_nodes_calls = 0;
    const size_t MERGE_CALL_BATCH_SIZE = 100000;
    double merge_nodes_call_time = 0.0;
    for (size_t t = 0; t < nb_threads; t++) {
        topology_workers.emplace_back([&]{
            vector<MergeCall> merge_calls;
            merge_calls.reserve(MERGE_CALL_BATCH_SIZE + 10);
            size_t merge_nodes_calls = 0;
            size_t offset = (dbg->size() / nb_threads) * t;
            size_t limit = (dbg->size() / nb_threads) * (t + 1);
            if (t == nb_threads - 1) {
                limit = dbg->size();
            }

            auto unitig_iterator = dbg->begin();

            for (size_t unitig_index = offset; unitig_index < limit; unitig_index++) {
                unitig_iterator.setIndex(unitig_index, h_kmer_ccov_orders);
                auto unitig = *unitig_iterator;

                for (const auto& successor: unitig.getSuccessors()) {
                    //cout << "successor.getIndex() = " << successor.getIndex(h_kmer_ccov_ranks) << endl;
                    MergeCall merge_call;
                    merge_call.u1_index = unitig_index;
                    merge_call.u1_strand = unitig.strand;
                    merge_call.u2_index = successor.getIndex(h_kmer_ccov_ranks);
                    merge_call.u2_strand = successor.strand;
                    merge_calls.push_back(merge_call);
                }
                for (const auto& predecessor: unitig.getPredecessors()) {
                    //cout << "predecessor.getIndex() = " << predecessor.getIndex(h_kmer_ccov_ranks) << endl;
                    MergeCall merge_call;
                    merge_call.u1_index = predecessor.getIndex(h_kmer_ccov_ranks);
                    merge_call.u1_strand = predecessor.strand;
                    merge_call.u2_index = unitig_index;
                    merge_call.u2_strand = unitig.strand;
                    merge_calls.push_back(merge_call);
                }

                if (merge_calls.size() > MERGE_CALL_BATCH_SIZE) {
                    unique_lock<mutex> lock(rust_helsitigs_mutex);
                    merge_nodes_call_time += do_merge_calls(merge_calls);
                    merge_nodes_calls += merge_calls.size();
                    merge_calls.clear();
                }
            }

            unique_lock<mutex> lock(rust_helsitigs_mutex);
            merge_nodes_call_time += do_merge_calls(merge_calls);
            merge_nodes_calls += merge_calls.size();
            merge_calls.clear();
            total_merge_nodes_calls += merge_nodes_calls;
        });
    }

    for (auto& t : topology_workers) t.join();

    /*size_t merge_nodes_calls = 0;
    for (const auto unitig : *dbg) {
        auto unitig_index = unitig.getIndex(h_kmer_ccov_ranks);
        //cout << "unitig.getIndex() = " << unitig.getIndex(h_kmer_ccov_ranks) << endl;
        for (const auto& successor: unitig.getSuccessors()) {
            //cout << "successor.getIndex() = " << successor.getIndex(h_kmer_ccov_ranks) << endl;
            helsitigs_merge_nodes(unitig_index, unitig.strand, successor.getIndex(h_kmer_ccov_ranks), successor.strand);
            merge_nodes_calls += 1;
        }
        for (const auto& predecessor: unitig.getPredecessors()) {
            //cout << "predecessor.getIndex() = " << predecessor.getIndex(h_kmer_ccov_ranks) << endl;
            helsitigs_merge_nodes(predecessor.getIndex(h_kmer_ccov_ranks), predecessor.strand, unitig_index, unitig.strand);
            merge_nodes_calls += 1;
        }

        // len is length of the mapping in kmers
        unitig_weights.push_back(unitig.len);
    }*/
    const auto stop_merge_nodes = std::chrono::high_resolution_clock::now();
    double duration_merge_nodes = ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop_merge_nodes - start_merge_nodes).count()) / 1e6;
    cout << "Took " << std::fixed << std::setprecision(3) << duration_merge_nodes << "s for " << total_merge_nodes_calls << " to send topology using " << nb_threads << " threads, of this " << merge_nodes_call_time << "s are from mutual exclusive calls to helsitigs_merge_nodes" << std::endl;

    helsitigs_build_graph(unitig_weights.data());
    const auto stop_send = std::chrono::high_resolution_clock::now();
    double duration_send = ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop_send - start_send).count()) / 1e6;
    cout << "Took " << std::fixed << std::setprecision(3) << duration_send << "s to build topology with rust_helsitigs" << std::endl;

    const auto start_compute = std::chrono::high_resolution_clock::now();
    vector<ptrdiff_t> tigs_edge_out(dbg->size() * 2, 0);
    vector<size_t> tigs_insert_out(dbg->size() * 2, 0);
    vector<size_t> tigs_out_limits(dbg->size(), 0);
    helsitigs_compute_tigs(tigs, nb_threads, k_, matching_file_prefix.c_str(), tigs_edge_out.data(), tigs_insert_out.data(), tigs_out_limits.data());
    const auto stop_compute = std::chrono::high_resolution_clock::now();
    const double duration_compute = ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop_compute - start_compute).count()) / 1e6;
    cout << "Took " << std::fixed << std::setprecision(3) << duration_compute << "s to compute tigs in rust_helsitigs" << std::endl;

    const auto start_receive = std::chrono::high_resolution_clock::now();
    hmap_min_unitigs = MinimizerIndex(dbg->hmap_min_unitigs.sz());

    size_t nb_tigs = 0;
    for (const auto limit : tigs_out_limits) {
        if (limit != 0) {
            nb_tigs ++;
        } else {
            break;
        }
    }

    cout << "CompactedDBG::convert_tigs(): Retrieving tigs from rust_helsitigs" << endl;
    size_t offset = 0;
    double duration_coloring = 0.0;
    auto unitig_iterator = dbg->begin();
    //auto debug_unitig_iterator = dbg->begin();
    for (const auto limit : tigs_out_limits) {
        if (limit == 0) {
            break;
        }

        string tig;
        size_t last_insert = 0;
        for (size_t i = offset; i < limit; i++) {
            auto edge = tigs_edge_out[i];
            auto insert = tigs_insert_out[i];

            if (insert == 0) {
                unitig_iterator.setIndex(abs(edge), h_kmer_ccov_orders);
                /*debug_unitig_iterator.setIndex(abs(edge));
                if (unitig_iterator != debug_unitig_iterator) {
                    cout << "CompactedDBG::convert_tigs(): Error: unitig iterators differ\nunitig_iterator:       " << unitig_iterator << "\ndebug_unitig_iterator: " << debug_unitig_iterator << endl;
                    exit(1);
                }*/
                const auto unitig_mapping = *unitig_iterator;

                auto sequence = unitig_mapping.mappedSequenceToString();
                if (edge < 0) {
                    sequence = reverse_complement(sequence);
                }

                if (tig.empty()) {
                    tig = sequence;
                } else {
                    tig.append(sequence.begin() + k_ - 1 - last_insert, sequence.end());
                    last_insert = 0;
                }
            } else {
                if (last_insert != 0) {
                    cerr << "Error: tig contains two consecutive dummy edges" << endl;
                    return false;
                }

                last_insert = insert;

                if (tig.empty()) {
                    cerr << "Error: first tig edge is dummy" << endl;
                    return false;
                }
            }
        }

        if (last_insert != 0) {
            cerr << "Error: last tig edge is dummy" << endl;
            return false;
        }

        auto tiglen = tig.length();
        auto tigid = (tiglen == k_) ? km_unitigs.size() : v_unitigs.size();
        //cout << "Adding tig of length " << tiglen << " with id " << tigid << endl;
        bool isAbundant = addUnitig(tig, tigid);
        if (isAbundant) {
            auto km_rep = Kmer(tig.c_str()).rep();
            auto it = h_kmers_ccov.find(km_rep);
            if (it == h_kmers_ccov.end()) {
                cout << "Insertion of " << km_rep.toString() << " failed" << endl;
            }
            tigid = it.h;
        }
        UnitigMap<U, G> um(tigid, 0, tiglen - k_ + 1, tiglen, (tiglen == k_) && !isAbundant, isAbundant, true, this);
        if (isAbundant) {
            auto head_key = h_kmers_ccov.find(tigid);
            if (head_key == h_kmers_ccov.end()) {
                cout << "tigid is wrong" << endl;
            }
            auto head = um.getMappedHead();
        }

        if (um.isEmpty) {
            cerr << "um is empty" << endl;
            exit(1);
        }

        // Add colors
        const auto start_coloring = std::chrono::high_resolution_clock::now();
        colorUnitig(dbg,
                    um,
                    h_kmer_ccov_orders,
                    tigs_edge_out.data() + offset,
                    tigs_edge_out.data() + limit,
                    tigs_insert_out.data() + offset,
                    tigs_insert_out.data() + limit,
                    nb_threads,
                    offset == 0, nb_tigs);
        const auto stop_coloring = std::chrono::high_resolution_clock::now();
        duration_coloring += ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop_coloring - start_coloring).count()) / 1e6;

        offset = limit;
    }

    const auto stop_receive = std::chrono::high_resolution_clock::now();
    const double duration_receive = ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop_receive - start_receive).count()) / 1e6;
    cout << "Took " << std::fixed << std::setprecision(3) << duration_receive << "s to receive tigs from rust_helsitigs" << std::endl;
    cout << "Took " << std::fixed << std::setprecision(3) << duration_coloring << "s to transfer colors while receiving tigs from rust_helsitigs" << std::endl;


    cout << "CompactedDBG::convert_tigs(): end" << endl;
    return true;
}

template<typename U, typename G>
void CompactedDBG<U, G>::colorUnitig(const CompactedDBG<U,G>* dbg,
                                     const UnitigMap<U,G>& um,
                                     const vector<size_t>& h_kmers_ccov_orders,
                                     const ptrdiff_t* tigs_edge_out_offset,
                                     const ptrdiff_t* tigs_edge_out_limit,
                                     const size_t* tigs_insert_out_offset,
                                     const size_t* tigs_insert_out_limit,
                                     const size_t nb_threads,
                                     const bool initialise_data,
                                     const size_t total_tigs) {
    //cout << "CompactedDBG::colorUnitig(): do nothing, since we are not colored" << endl;
}

#endif