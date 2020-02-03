#ifndef BIFROST_COLOREDCDBG_TCC
#define BIFROST_COLOREDCDBG_TCC

template<typename U>
ColoredCDBG<U>::ColoredCDBG(int kmer_length, int minimizer_length) : CompactedDBG<DataAccessor<U>, DataStorage<U>>(kmer_length, minimizer_length){

    invalid = this->isInvalid();
}

template<typename U>
ColoredCDBG<U>::ColoredCDBG(const ColoredCDBG& o) : CompactedDBG<DataAccessor<U>, DataStorage<U>>(o), invalid(o.invalid) {}

template<typename U>
ColoredCDBG<U>::ColoredCDBG(ColoredCDBG&& o) :  CompactedDBG<DataAccessor<U>, DataStorage<U>>(move(o)), invalid(o.invalid) {}

template<typename U>
void ColoredCDBG<U>::clear(){

    invalid = true;

    this->getData()->clear();
    CompactedDBG<DataAccessor<U>, DataStorage<U>>::clear();
}

template<typename U>
ColoredCDBG<U>& ColoredCDBG<U>::operator=(const ColoredCDBG& o) {

    CompactedDBG<DataAccessor<U>, DataStorage<U>>::operator=(o);

    invalid = o.invalid;

    return *this;
}

template<typename U>
ColoredCDBG<U>& ColoredCDBG<U>::operator=(ColoredCDBG&& o) {

    if (this != &o) {

        CompactedDBG<DataAccessor<U>, DataStorage<U>>::operator=(move(o));

        invalid = o.invalid;

        o.clear();
    }

    return *this;
}

template<typename U>
ColoredCDBG<U>& ColoredCDBG<U>::operator+=(const ColoredCDBG& o) {

    if (this != &o) merge(o, 1, false);

    return *this;
}

template<typename U>
bool ColoredCDBG<U>::operator==(const ColoredCDBG& o) const {

    if (!invalid && !this->isInvalid() && !o.isInvalid() && (this->getK() == o.getK()) && (this->size() == o.size())){

        const size_t k = this->getK();

        for (const auto& unitig : *this){

            const_UnitigColorMap<U> unitig_o(o.find(unitig.getUnitigHead(), true));

            if (unitig_o.isEmpty) return false;
            else {

                unitig_o.dist = 0;
                unitig_o.len = unitig_o.size - k + 1;

                const string unitig_o_str = unitig_o.strand ? unitig_o.referenceUnitigToString() : reverse_complement(unitig_o.referenceUnitigToString());

                if (unitig_o_str != unitig.referenceUnitigToString()) return false;
                else {

                    const UnitigColors* uc = unitig.getData()->getUnitigColors(unitig);
                    const UnitigColors* uc_o = unitig_o.getData()->getUnitigColors(unitig_o);

                    if ((uc != nullptr) && (uc_o != nullptr)){

                        if (!uc->isEqual(unitig, *uc_o, unitig_o)) return false;
                    }
                    else if ((uc != nullptr) != (uc_o != nullptr)) return false;
                }
            }
        }

        return true;
    }

    return false;
}

template<typename U>
inline bool ColoredCDBG<U>::operator!=(const ColoredCDBG& o) const {

    return !operator==(o);
}

template<typename U>
bool ColoredCDBG<U>::merge(const ColoredCDBG& o, const size_t nb_threads, const bool verbose){

    bool ret = true;

    if (invalid){

         if (verbose) cerr << "ColoredCDBG::merge(): Current graph is invalid." << endl;
         ret = false;
    }

    if (o.invalid){

         if (verbose) cerr << "ColoredCDBG::merge(): Graph to merge is invalid." << endl;
         ret = false;
    }

    if (this->getK() != o.getK()){

         if (verbose) cerr << "ColoredCDBG::merge(): The graphs to merge do not have the same k-mer length." << endl;
         ret = false;
    }

    if (this == &o){

         if (verbose) cerr << "ColoredCDBG::merge(): Cannot merge graph with itself." << endl;
         ret = false;
    }

    if (ret){

        const size_t sz_before = this->size();

        for (auto& unitig : *this) unitig.setFullCoverage();

        ret = CompactedDBG<DataAccessor<U>, DataStorage<U>>::annotateSplitUnitigs(o, nb_threads, verbose);

        if (ret){

            const size_t sz_after = this->size();
            const pair<size_t, size_t> p1 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::getSplitInfoAllUnitigs();

            resizeDataUC(sz_after + (p1.second - p1.first), nb_threads);

            const pair<size_t, size_t> p2 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::splitAllUnitigs();
            const size_t joined = (p1.second != 0) ? CompactedDBG<DataAccessor<U>, DataStorage<U>>::joinUnitigs() : 0;

            if (verbose){

                cout << "CompactedDBG::merge(): Added " << (sz_after - sz_before) << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Split " << p1.first << " unitigs into " << p1.second << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Joined " << joined << " unitigs." << endl;
                cout << "CompactedDBG::merge(): " << this->size() << " unitigs after merging." << endl;
            }

            for (size_t i = 0; i < o.getNbColors(); ++i) this->getData()->color_names.push_back(o.getColorName(i));

            return CompactedDBG<DataAccessor<U>, DataStorage<U>>::mergeData(o, nb_threads, verbose);
        }
    }

    return false;
}

template<typename U>
bool ColoredCDBG<U>::merge(ColoredCDBG&& o, const size_t nb_threads, const bool verbose){

    bool ret = true;

    if (invalid){

         if (verbose) cerr << "ColoredCDBG::merge(): Current graph is invalid." << endl;
         ret = false;
    }

    if (o.invalid){

         if (verbose) cerr << "ColoredCDBG::merge(): Graph to merge is invalid." << endl;
         ret = false;
    }

    if (this->getK() != o.getK()){

         if (verbose) cerr << "ColoredCDBG::merge(): The graphs to merge do not have the same k-mer length." << endl;
         ret = false;
    }

    if (this == &o){

         if (verbose) cerr << "ColoredCDBG::merge(): Cannot merge graph with itself." << endl;
         ret = false;
    }

    if (ret){

        const size_t sz_before = this->size();

        for (auto& unitig : *this) unitig.setFullCoverage();

        ret = CompactedDBG<DataAccessor<U>, DataStorage<U>>::annotateSplitUnitigs(o, nb_threads, verbose);

        if (ret){

            const size_t sz_after = this->size();
            const pair<size_t, size_t> p1 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::getSplitInfoAllUnitigs();

            resizeDataUC(sz_after + (p1.second - p1.first), nb_threads);

            const pair<size_t, size_t> p2 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::splitAllUnitigs();
            const size_t joined = (p1.second != 0) ? CompactedDBG<DataAccessor<U>, DataStorage<U>>::joinUnitigs() : 0;

            if (verbose){

                cout << "CompactedDBG::merge(): Added " << (sz_after - sz_before) << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Split " << p1.first << " unitigs into " << p1.second << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Joined " << joined << " unitigs." << endl;
                cout << "CompactedDBG::merge(): " << this->size() << " unitigs after merging." << endl;
            }

            for (size_t i = 0; i < o.getNbColors(); ++i) this->getData()->color_names.push_back(o.getColorName(i));

            const bool ret = CompactedDBG<DataAccessor<U>, DataStorage<U>>::mergeData(move(o), nb_threads, verbose);

            o.clear();

            return ret;
        }
    }

    return false;
}

template<typename U>
bool ColoredCDBG<U>::merge(const vector<ColoredCDBG>& v, const size_t nb_threads, const bool verbose){

    bool ret = true;

    if (invalid){

         if (verbose) cerr << "ColoredCDBG::merge(): Current graph is invalid." << endl;
         ret = false;
    }

    for (const auto& ccdbg : v){

        if (ccdbg.invalid){

             if (verbose) cerr << "ColoredCDBG::merge(): One of the graph to merge is invalid." << endl;
             ret = false;
        }

        if (this->getK() != ccdbg.getK()){

             if (verbose) cerr << "ColoredCDBG::merge(): The graphs to merge do not have the same k-mer length." << endl;
             ret = false;
        }

        if (this == &ccdbg){

             if (verbose) cerr << "ColoredCDBG::merge(): Cannot merge graph with itself." << endl;
             ret = false;
        }
    }

    if (ret){

        const size_t sz_before = this->size();

        for (auto& unitig : *this) unitig.setFullCoverage();

        for (const auto& ccdbg : v){

            ret = CompactedDBG<DataAccessor<U>, DataStorage<U>>::annotateSplitUnitigs(ccdbg, nb_threads, verbose);

            if (!ret) break;
        }

        if (ret){

            const size_t sz_after = this->size();
            const pair<size_t, size_t> p1 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::getSplitInfoAllUnitigs();

            resizeDataUC(sz_after + (p1.second - p1.first), nb_threads);

            const pair<size_t, size_t> p2 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::splitAllUnitigs();
            const size_t joined = (p1.second != 0) ? CompactedDBG<DataAccessor<U>, DataStorage<U>>::joinUnitigs() : 0;

            if (verbose){

                cout << "CompactedDBG::merge(): Added " << (sz_after - sz_before) << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Split " << p1.first << " unitigs into " << p1.second << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Joined " << joined << " unitigs." << endl;
                cout << "CompactedDBG::merge(): " << this->size() << " unitigs after merging." << endl;
            }

            for (const auto& ccdbg : v){

                for (size_t i = 0; i < ccdbg.getNbColors(); ++i) this->getData()->color_names.push_back(ccdbg.getColorName(i));

                if (!CompactedDBG<DataAccessor<U>, DataStorage<U>>::mergeData(ccdbg, nb_threads, verbose)) return false;
            }

            return true;
        }
    }

    return false;
}

template<typename U>
bool ColoredCDBG<U>::merge(vector<ColoredCDBG>&& v, const size_t nb_threads, const bool verbose){

    bool ret = true;

    if (invalid){

         if (verbose) cerr << "ColoredCDBG::merge(): Current graph is invalid." << endl;
         ret = false;
    }

    for (const auto& ccdbg : v){

        if (ccdbg.invalid){

             if (verbose) cerr << "ColoredCDBG::merge(): One of the graph to merge is invalid." << endl;
             ret = false;
        }

        if (this->getK() != ccdbg.getK()){

             if (verbose) cerr << "ColoredCDBG::merge(): The graphs to merge do not have the same k-mer length." << endl;
             ret = false;
        }

        if (this == &ccdbg){

             if (verbose) cerr << "ColoredCDBG::merge(): Cannot merge graph with itself." << endl;
             ret = false;
        }
    }

    if (ret){

        const size_t sz_before = this->size();

        for (auto& unitig : *this) unitig.setFullCoverage();

        for (const auto& ccdbg : v){

            ret = CompactedDBG<DataAccessor<U>, DataStorage<U>>::annotateSplitUnitigs(ccdbg, nb_threads, verbose);

            if (!ret) break;
        }

        if (ret){

            const size_t sz_after = this->size();
            const pair<size_t, size_t> p1 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::getSplitInfoAllUnitigs();

            resizeDataUC(sz_after + (p1.second - p1.first), nb_threads);

            const pair<size_t, size_t> p2 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::splitAllUnitigs();
            const size_t joined = (p1.second != 0) ? CompactedDBG<DataAccessor<U>, DataStorage<U>>::joinUnitigs() : 0;

            if (verbose){

                cout << "CompactedDBG::merge(): Added " << (sz_after - sz_before) << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Split " << p1.first << " unitigs into " << p1.second << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Joined " << joined << " unitigs." << endl;
                cout << "CompactedDBG::merge(): " << this->size() << " unitigs after merging." << endl;
            }

            for (auto& ccdbg : v){

                for (size_t i = 0; i < ccdbg.getNbColors(); ++i) this->getData()->color_names.push_back(ccdbg.getColorName(i));

                if (!CompactedDBG<DataAccessor<U>, DataStorage<U>>::mergeData(move(ccdbg), nb_threads, verbose)) return false;

                ccdbg.clear();
            }

            return true;
        }
    }

    return false;
}

template<typename U>
bool ColoredCDBG<U>::buildGraph(const CCDBG_Build_opt& opt){

    if (!invalid){

        CDBG_Build_opt opt_ = opt;

        invalid = !this->build(opt_);
    }
    else cerr << "ColoredCDBG::buildGraph(): Graph is invalid and cannot be built." << endl;

    return !invalid;
}

template<typename U>
bool ColoredCDBG<U>::buildColors(const CCDBG_Build_opt& opt){

    if (!invalid){

        initUnitigColors(opt);
        buildUnitigColors(opt.nb_threads);
    }
    else cerr << "ColoredCDBG::buildColors(): Graph is invalid (maybe not built yet?) and colors cannot be mapped." << endl;

    return !invalid;
}

template<typename U>
bool ColoredCDBG<U>::write(const string& prefix_output_filename, const size_t nb_threads, const bool verbose) const {

    if (!CompactedDBG<DataAccessor<U>, DataStorage<U>>::write(prefix_output_filename, nb_threads, true, verbose)) return false; // Write graph

    return this->getData()->write(prefix_output_filename, verbose); // Write colors
}

template<typename U>
bool ColoredCDBG<U>::read(const string& input_graph_filename, const string& input_colors_filename, const size_t nb_threads, const bool verbose) {

    bool valid_input_files = true;

    if (input_graph_filename.length() != 0){

        if (check_file_exists(input_graph_filename)){

            FILE* fp = fopen(input_graph_filename.c_str(), "r");

            if (fp == NULL) {

                cerr << "ColoredCDBG::read(): Could not open input graph file " << input_graph_filename << endl;
                valid_input_files = false;
            }
            else fclose(fp);
        }
        else {

            cerr << "ColoredCDBG::read(): Input graph file " << input_graph_filename << " does not exist." << endl;
            valid_input_files = false;
        }
    }
    else {

        cerr << "ColoredCDBG::read(): No input graph file provided." << endl;
        valid_input_files = false;
    }

    if (input_colors_filename.length() != 0){

        if (check_file_exists(input_colors_filename)){

            FILE* fp = fopen(input_colors_filename.c_str(), "rb");

            if (fp == NULL) {

                cerr << "ColoredCDBG::read(): Could not open input colors file " << input_colors_filename << endl;
                valid_input_files = false;
            }
            else fclose(fp);
        }
        else {

            cerr << "ColoredCDBG::read(): Input colors file " << input_colors_filename << " does not exist." << endl;
            valid_input_files = false;
        }
    }
    else {

        cerr << "ColoredCDBG::read(): No input colors file provided." << endl;
        valid_input_files = false;
    }

    if (valid_input_files){

        if (verbose) cout << "ColoredCDBG::read(): Reading graph." << endl;

        if (!CompactedDBG<DataAccessor<U>, DataStorage<U>>::read(input_graph_filename, nb_threads, verbose)) return false; // Read graph

        if (verbose) cout << "ColoredCDBG::read(): Reading colors." << endl;

        if (!this->getData()->read(input_colors_filename, nb_threads, verbose)) return false; // Read colors

        if (verbose) cout << "ColoredCDBG::read(): Joining unitigs to their color sets." << endl;

        GFA_Parser graph(input_graph_filename);

        graph.open_read();

        auto reading_function = [&graph](vector<pair<Kmer, uint8_t>>& unitig_tags, const size_t chunk_size) {

            size_t i = 0;
            size_t graph_file_id = 0;

            bool new_file_opened = false;

            GFA_Parser::GFA_line r = graph.read(graph_file_id, new_file_opened, true);

            while ((i < chunk_size) && ((r.first != nullptr) || (r.second != nullptr))){

                if (r.first != nullptr){ // It is a sequence

                    if (r.first->tags.empty()){

                        cerr << "ColoredCDBG::read(): One sequence line in GFA file has no DataAccessor tag. Operation aborted." << endl;
                        return false;
                    }

                    size_t i = 0;

                    for (; i < r.first->tags.size(); ++i){

                        if (r.first->tags[i].substr(0, 5) == "DA:Z:") break;
                    }

                    if (i == r.first->tags.size()){

                        cerr << "ColoredCDBG::read(): One sequence line in GFA file has no DataAccessor tag. Operation aborted." << endl;
                        return false;
                    }

                    unitig_tags.push_back({Kmer(r.first->seq.c_str()), atoi(r.first->tags[i].c_str() + 5)});

                    ++i;
                }

                r = graph.read(graph_file_id, new_file_opened, true);
            }

            return ((r.first != nullptr) || (r.second != nullptr));
        };

        auto join_function = [this](const vector<pair<Kmer, uint8_t>>& unitig_tags) {

            for (const auto& p : unitig_tags){

                UnitigColorMap<U> ucm(this->find(p.first, true));

                if (ucm.isEmpty){

                    cerr << "ColoredCDBG::read(): Internal error, operation aborted." << endl;
                    cerr << "ColoredCDBG::read(): A unitig from GFA file is not found in the in-memory graph." << endl;
                    cerr << "ColoredCDBG::read(): Graph from GFA file possibly incorrectly compacted." << endl;

                    return false;
                }

                DataAccessor<U>* da = ucm.getData();

                *da = DataAccessor<U>(p.second);

                if (!ucm.strand){ // Unitig has been inserted in reverse-complement, need to reverse order of color sets

                    UnitigColors* uc = da->getUnitigColors(ucm);

                    UnitigColors r_uc = uc->reverse(ucm);

                    *uc = move(r_uc);
                }
            }

            return true;
        };

        {
            const size_t chunk = 10000;

            vector<thread> workers; // need to keep track of threads so we can join them
            vector<vector<pair<Kmer, uint8_t>>> v(nb_threads);

            mutex mutex_file;

            bool file_valid_for_read = true;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        while (true) {

                            {
                                unique_lock<mutex> lock(mutex_file);

                                if (!file_valid_for_read) return;

                                file_valid_for_read = reading_function(v[t], chunk);

                            }

                            join_function(v[t]);
                            v[t].clear();
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }
    }

    return valid_input_files;
}

template<typename U>
void ColoredCDBG<U>::initUnitigColors(const CCDBG_Build_opt& opt, const size_t max_nb_hash){

    vector<string> v_files(opt.filename_seq_in);

    v_files.insert(v_files.end(), opt.filename_ref_in.begin(), opt.filename_ref_in.end());

    DataStorage<U>* ds = this->getData();
    DataStorage<U> new_ds(max_nb_hash, this->size(), v_files);

    *ds = move(new_ds);

    v_files.clear();

    const size_t chunk = 1000;

    vector<thread> workers; // need to keep track of threads so we can join them

    typename ColoredCDBG<U>::iterator g_a = this->begin();
    typename ColoredCDBG<U>::iterator g_b = this->end();

    mutex mutex_it;

    for (size_t t = 0; t < opt.nb_threads; ++t){

        workers.emplace_back(

            [&, t]{

                typename ColoredCDBG<U>::iterator l_a, l_b;

                while (true) {

                    {
                        unique_lock<mutex> lock(mutex_it);

                        if (g_a == g_b) return;

                        l_a = g_a;
                        l_b = g_a;

                        for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                        g_a = l_b;
                    }

                    for (auto& it_unitig = l_a; it_unitig != l_b; ++it_unitig) {

                        *(it_unitig->getData()) = ds->insert(*it_unitig).first;
                    }
                }
            }
        );
    }

    for (auto& t : workers) t.join();

    //cout << "Number of unitigs not hashed is " << ds->overflow.size() << " on " << ds->nb_cs << " unitigs." << endl;
}

template<typename U>
void ColoredCDBG<U>::resizeDataUC(const size_t sz, const size_t nb_threads, const size_t max_nb_hash){

    DataStorage<U>* ds = this->getData();

    DataStorage<U> new_ds(max_nb_hash, sz, ds->color_names);

    const size_t chunk = 100;

    vector<thread> workers; // need to keep track of threads so we can join them

    typename ColoredCDBG<U>::iterator g_a = this->begin();
    typename ColoredCDBG<U>::iterator g_b = this->end();

    mutex mutex_it;

    for (size_t t = 0; t < nb_threads; ++t){

        workers.emplace_back(

            [&, t]{

                typename ColoredCDBG<U>::iterator l_a, l_b;

                while (true) {

                    {
                        unique_lock<mutex> lock(mutex_it);

                        if (g_a == g_b) return;

                        l_a = g_a;
                        l_b = g_a;

                        for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                        g_a = l_b;
                    }

                    for (auto& it_unitig = l_a; it_unitig != l_b; ++it_unitig) {

                        UnitigColors* uc = ds->getUnitigColors(*it_unitig);
                        U* data = ds->getData(*it_unitig);

                        if ((uc != nullptr) || (data != nullptr)){

                            const pair<DataAccessor<U>, pair<UnitigColors*, U*>> p  = new_ds.insert(*it_unitig);

                            *(it_unitig->getData()) = p.first;

                            if (uc != nullptr) *(p.second.first) = move(*uc);
                            if (data != nullptr) *(p.second.second) = move(*data);
                        }
                    }
                }
            }
        );
    }

    for (auto& t : workers) t.join();

    *ds = move(new_ds);

    //cout << "Number of unitigs not hashed is " << ds->overflow.size() << " on " << ds->nb_cs << " unitigs." << endl;
}

template<>
inline void ColoredCDBG<void>::resizeDataUC(const size_t sz, const size_t nb_threads, const size_t max_nb_hash){

    DataStorage<void>* ds = this->getData();

    DataStorage<void> new_ds(max_nb_hash, sz, ds->color_names);

    const size_t chunk = 100;

    vector<thread> workers; // need to keep track of threads so we can join them

    typename ColoredCDBG<void>::iterator g_a = this->begin();
    typename ColoredCDBG<void>::iterator g_b = this->end();

    mutex mutex_it;

    for (size_t t = 0; t < nb_threads; ++t){

        workers.emplace_back(

            [&, t]{

                typename ColoredCDBG<void>::iterator l_a, l_b;

                while (true) {

                    {
                        unique_lock<mutex> lock(mutex_it);

                        if (g_a == g_b) return;

                        l_a = g_a;
                        l_b = g_a;

                        for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                        g_a = l_b;
                    }

                    for (auto& it_unitig = l_a; it_unitig != l_b; ++it_unitig) {

                        UnitigColors* uc = ds->getUnitigColors(*it_unitig);

                        if (uc != nullptr){

                            const pair<DataAccessor<void>, pair<UnitigColors*, void*>> p  = new_ds.insert(*it_unitig);

                            *(it_unitig->getData()) = p.first;
                            *(p.second.first) = move(*uc);
                        }
                    }
                }
            }
        );
    }

    for (auto& t : workers) t.join();

    *ds = move(new_ds);

    //cout << "Number of unitigs not hashed is " << ds->overflow.size() << " on " << ds->nb_cs << " unitigs." << endl;
}

template<typename U>
void ColoredCDBG<U>::buildUnitigColors(const size_t nb_threads){

    DataStorage<U>* ds = this->getData();

    const int k_ = this->getK();

    const size_t nb_locks = nb_threads * 1024;
    const size_t chunk_size = 64;
    const size_t max_len_seq = 1024;
    const size_t thread_seq_buf_sz = chunk_size * max_len_seq;
    const size_t thread_col_buf_sz = (thread_seq_buf_sz / (k_ + 1)) + 1;

    size_t prev_file_id = 0;

    size_t pos_read = 0;
    size_t len_read = 0;

    bool next_file = true;

    string s;

    FileParser fp(ds->color_names);

    std::atomic_flag* cs_locks = new std::atomic_flag[nb_locks];

    for (size_t i = 0; i < nb_locks; ++i) cs_locks[i].clear();

    // Main worker thread
    auto worker_function = [&](char* seq_buf, const size_t seq_buf_sz, const size_t* col_buf) {

        char* str = seq_buf;
        const char* str_end = &seq_buf[seq_buf_sz];

        size_t c_id = 0;

        while (str < str_end) { // for each input

            const int len = strlen(str);

            for (char* s = str; s != &str[len]; ++s) *s &= 0xDF;

            for (size_t i = 0; i < len - k_ + 1; i += max_len_seq - k_ + 1){

                const int curr_len = min(len - i, max_len_seq);
                const char saved_char = str[i + curr_len];
                const char* str_tmp = &str[i];

                str[i + curr_len] = '\0';

                for (KmerIterator it_km(str_tmp), it_km_end; it_km != it_km_end; ++it_km) {

                    UnitigColorMap<U> um = this->find(it_km->first);

                    if (!um.isEmpty) {

                        if (um.strand || (um.dist != 0)){

                            um.len = 1 + um.lcp(str_tmp, it_km->second + k_, um.strand ? um.dist + k_ : um.dist - 1, !um.strand);

                            //if ((um.size != k_) && !um.strand) um.dist -= um.len - 1;
                            um.dist -= (um.len - 1) & (static_cast<size_t>((um.size == k_) || um.strand) - 1);

                            it_km += um.len - 1;
                        }

                        const uint64_t id_lock = ds->getHash(um) % nb_locks;
                        UnitigColors* uc = ds->getUnitigColors(um);

                        while (cs_locks[id_lock].test_and_set(std::memory_order_acquire)); // Set the corresponding lock

                        uc->add(um, col_buf[c_id]);

                        cs_locks[id_lock].clear(std::memory_order_release);
                    }
                }

                str[i + curr_len] = saved_char;
            }

            str += len + 1;
            ++c_id;
        }
    };

    auto reading_function = [&](char* seq_buf, size_t& seq_buf_sz, size_t* col_buf) {

        size_t file_id = prev_file_id;
        size_t i = 0;

        const size_t sz_buf = thread_seq_buf_sz - k_;

        const char* s_str = s.c_str();

        seq_buf_sz = 0;

        while (seq_buf_sz < sz_buf) {

            const bool new_reading = (pos_read >= len_read);

            if (!new_reading || fp.read(s, file_id)) {

                //pos_read = (new_reading ? 0 : pos_read);
                pos_read &= static_cast<size_t>(new_reading) - 1;

                len_read = s.length();
                s_str = s.c_str();

                if (len_read >= k_){

                    if ((thread_seq_buf_sz - seq_buf_sz - 1) < (len_read - pos_read)){

                        strncpy(&seq_buf[seq_buf_sz], &s_str[pos_read], thread_seq_buf_sz - seq_buf_sz - 1);

                        seq_buf[thread_seq_buf_sz - 1] = '\0';
                        col_buf[i++] = file_id;

                        pos_read += sz_buf - seq_buf_sz;
                        seq_buf_sz = thread_seq_buf_sz;

                        break;
                    }
                    else {

                        strcpy(&seq_buf[seq_buf_sz], &s_str[pos_read]);

                        col_buf[i++] = file_id;

                        seq_buf_sz += (len_read - pos_read) + 1;
                        pos_read = len_read;
                    }
                }
                else pos_read = len_read;
            }
            else {

                next_file = false;

                return true;
            }
        }

        const bool ret = (file_id != prev_file_id);

        next_file = true;
        prev_file_id = file_id;

        return ret;
    };

    {
        bool stop = false;

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_file;

        size_t prev_uc_sz = getCurrentRSS();

        char* buffer_seq = new char[nb_threads * thread_seq_buf_sz];
        size_t* buffer_seq_sz = new size_t[nb_threads];
        size_t* buffer_col = new size_t[nb_threads * thread_col_buf_sz];

        while (next_file){

            stop = false;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        while (true) {

                            {
                                unique_lock<mutex> lock(mutex_file);

                                if (stop) return;

                                stop = reading_function(&buffer_seq[t * thread_seq_buf_sz], buffer_seq_sz[t], &buffer_col[t * thread_col_buf_sz]);
                            }

                            worker_function(&buffer_seq[t * thread_seq_buf_sz], buffer_seq_sz[t], &buffer_col[t * thread_col_buf_sz]);
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();

            workers.clear();

            const size_t curr_uc_sz = getCurrentRSS();

            if ((curr_uc_sz - prev_uc_sz) >= 1073741824ULL){

                const size_t chunk = 1000;

                typename ColoredCDBG<U>::iterator g_a = this->begin();
                typename ColoredCDBG<U>::iterator g_b = this->end();

                mutex mutex_it;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&, t]{

                            typename ColoredCDBG<U>::iterator l_a, l_b;

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_it);

                                    if (g_a == g_b) return;

                                    l_a = g_a;
                                    l_b = g_a;

                                    for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                                    g_a = l_b;
                                }

                                while (l_a != l_b){

                                    l_a->getData()->getUnitigColors(*l_a)->optimizeFullColors(*l_a);
                                    ++l_a;
                                }
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();

                workers.clear();

                prev_uc_sz = getCurrentRSS();
            }
        }

        delete[] buffer_seq;
        delete[] buffer_seq_sz;
        delete[] buffer_col;
    }

    fp.close();

    //checkColors(ds->color_names);

    /*typedef std::unordered_map<uint64_t, pair<int64_t, size_t>> uc_unordered_map;

    uc_unordered_map u_map;

    mutex mutex_u_map;

    vector<thread> workers;

    auto add_hash_function = [&](typename ColoredCDBG<U>::iterator it_a, typename ColoredCDBG<U>::iterator it_b) {

        while (it_a != it_b) {

            const const_UnitigColorMap<U> unitig(*it_a);
            const UnitigColors* uc = unitig.getData()->getUnitigColors(unitig);
            const UnitigColors uc_full = uc->makeFullColors(unitig);
            const UnitigColors* uc_full_array = uc_full.getFullColorsPtr();

            if (uc_full_array[0].size() != 0){

                const pair<int64_t, size_t> pv(0 - static_cast<int64_t>(uc_full_array[0].getSizeInBytes()) - static_cast<int64_t>(sizeof(size_t)), 0);

                const int64_t to_add = static_cast<int64_t>(uc->getSizeInBytes());
                const int64_t to_rm = (static_cast<int64_t>(uc_full_array[1].getSizeInBytes() + 2 * sizeof(UnitigColors)));

                {
                    unique_lock<mutex> lock(mutex_u_map);

                    pair<uc_unordered_map::iterator, bool> p = u_map.insert(make_pair(uc_full_array[0].hash(), pv));

                    p.first->second.first += to_add - to_rm;
                }
            }

            ++it_a;
        }
    };

    auto add_shared_function = [&](typename ColoredCDBG<U>::iterator it_a, typename ColoredCDBG<U>::iterator it_b) {

        while (it_a != it_b) {

            const UnitigColorMap<U> unitig(*it_a);

            UnitigColors* uc = unitig.getData()->getUnitigColors(unitig);
            UnitigColors uc_full = uc->makeFullColors(unitig);
            UnitigColors* uc_full_array = uc_full.getFullColorsPtr();

            if (uc_full_array[0].size() != 0){

                uc_unordered_map::const_iterator it = u_map.find(uc_full_array[0].hash());

                if (it->second.first > 0){ // If there is some sharing

                    const size_t id_shared = it->second.second;
                    const uint64_t id_lock = id_shared % nb_locks;

                    bool move_full = false;

                    while (cs_locks[id_lock].test_and_set(std::memory_order_acquire)); // Set the corresponding lock

                    if (ds->shared_color_sets[id_shared].second == 0){

                        ds->shared_color_sets[id_shared].first = move(uc_full_array[0]);
                        ds->shared_color_sets[id_shared].second = 0;

                        uc_full_array[0] = ds->shared_color_sets[id_shared];
                        move_full = true;
                    }
                    else if (uc_full_array[0] == ds->shared_color_sets[id_shared]) {

                        uc_full_array[0] = ds->shared_color_sets[id_shared];
                        move_full = true;
                    }

                    cs_locks[id_lock].clear(std::memory_order_release);

                    if (move_full) *uc = move(uc_full);
                }
            }

            ++it_a;
        }
    };

    {
        const size_t chunk = 1000;

        typename ColoredCDBG<U>::iterator g_a = this->begin();
        typename ColoredCDBG<U>::iterator g_b = this->end();

        mutex mutex_it;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    typename ColoredCDBG<U>::iterator l_a, l_b;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_it);

                            if (g_a == g_b) return;

                            l_a = g_a;
                            l_b = g_a;

                            for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                            g_a = l_b;
                        }

                        add_hash_function(l_a, l_b);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        workers.clear();
    }

    for (auto& p : u_map){

        if (p.second.first > 0){

            p.second.second = ds->sz_shared_cs;
            ++(ds->sz_shared_cs);
        }
    }

    ds->shared_color_sets = new UnitigColors::SharedUnitigColors[ds->sz_shared_cs];

    {
        const size_t chunk = 1000;

        typename ColoredCDBG<U>::iterator g_a = this->begin();
        typename ColoredCDBG<U>::iterator g_b = this->end();

        mutex mutex_it;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    typename ColoredCDBG<U>::iterator l_a, l_b;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_it);

                            if (g_a == g_b) return;

                            l_a = g_a;
                            l_b = g_a;

                            for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                            g_a = l_b;
                        }

                        add_shared_function(l_a, l_b);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        workers.clear();
    }*/

    delete[] cs_locks;
}

template<typename U>
string ColoredCDBG<U>::getColorName(const size_t color_id) const {

    if (invalid){

        cerr << "ColoredCDBG::getColorName(): Graph is invalid or colors are not yet mapped to unitigs." << endl;
        return string();
    }

    const DataStorage<U>* ds = this->getData();

    if (color_id >= ds->color_names.size()){

        cerr << "ColoredCDBG::getColorName(): Color ID " << color_id << " is invalid, graph only has " <<
        ds->color_names.size() << " colors." << endl;

        return string();
    }

    return ds->color_names[color_id];
}

template<typename U>
vector<string> ColoredCDBG<U>::getColorNames() const {

    if (invalid){

        cerr << "ColoredCDBG::getColorNames(): Graph is invalid or colors are not yet mapped to unitigs." << endl;
        return vector<string>();
    }

    return this->getData()->color_names;
}


/*template<typename U>
void ColoredCDBG<U>::searchSequence(const string& seq, const bool exact, const bool insertion, const bool deletion, const bool substitution) const {

    if (invalid){

        cerr << "ColoredCDBG::searchSequence(): Graph is invalid and cannot be searched" << endl;

        return;
    }

    if (seq.length() < k_){

        cerr << "ColoredCDBG::searchSequence(): Query length is shorter than k-mer size" << endl;

        return;
    }

    const size_t nb_colors = getNbColors();

    string seqs;

    vector<Roaring> color_occ_r(nb_colors);
    vector<uint32_t> color_occ_u(nb_colors);

    auto worker_func = [&](const bool subst, const bool ins, const bool del, const size_t shift){

        const bool subst_or_ind = (subst || ins);
        const bool inexact = (subst_or_ind || del);

        const size_t ins_mask = static_cast<size_t>(!ins) - 1;
        const size_t del_mask = static_cast<size_t>(!del) - 1;

        const size_t end = 1ULL << (static_cast<size_t>(subst_or_ind) << 2);
        const size_t seq_len = seq.length();

        auto processUnitigMap = [&](const_UnitigColorMap<U> um, const size_t pos_seq){

            if (um.strand){

                size_t j = um.dist;

                const size_t end_pos = um.dist + um.len;

                while (j < end_pos){

                    size_t l_pos_seq = pos_seq + j - um.dist;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if (l_pos_seq + k_ - 1 >= seq_len) break;
                    else ++j;
                }

                um.len = j - um.dist;

                const UnitigColors* uc = um.getData()->getUnitigColors(um);

                UnitigColors::const_iterator it_uc = uc->begin(um);
                UnitigColors::const_iterator it_uc_end = uc->end();

                if (inexact){

                    for (; it_uc != it_uc_end; ++it_uc) color_occ_r[it_uc.getColorID()].add(it_uc.getKmerPosition() - um.dist + p.first);
                }
                else {

                    for (; it_uc != it_uc_end; ++it_uc) color_occ_u[it_uc.getColorID()] += 1;
                }
            }
            else {

                const size_t end_pos = um.dist + um.len;

                size_t j = end_pos - 1;

                while (j >= um.dist){

                    size_t l_pos_seq = pos_seq + um.dist + um.len - j - 1;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if (l_pos_seq + k_ - 1 >= seq_len) break;
                    else --j;
                }

                um.dist = j + 1;
                um.len = end_pos - j;

                const UnitigColors* uc = um.getData()->getUnitigColors(um);

                UnitigColors::const_iterator it_uc = uc->begin(um);
                UnitigColors::const_iterator it_uc_end = uc->end();

                if (inexact){

                    for (; it_uc != it_uc_end; ++it_uc) color_occ_r[it_uc.getColorID()].add(max_pos_um - it_uc.getKmerPosition() + p.first);
                }
                else {

                    for (; it_uc != it_uc_end; ++it_uc) color_occ_u[it_uc.getColorID()] += 1;
                }
            }
        };

        for (size_t i = 0; i != end; ++i){

            if (subst_or_ind){

                for (size_t j = shift; j < seqs.length(); j += k_) seqs[j] = alpha[i];
            }

            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(seqs.c_str(), seqs.length(), k_, g_, RepHash(), true);
            minHashResultIterator<RepHash> it_min = *mhi, it_min_end;
            minHashResult mhr = *it_min;

            Minimizer minz = Minimizer(seqs.c_str() + mhr.pos).rep();

            pair<size_t, bool> minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

            for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                mhr = *it_min;
                minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
            }

            size_t pos_seq = 0;
            size_t l_pos_seq = 0;
            size_t shift_pos_seq = 0;

            shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

            l_pos_seq -= (ins_mask & shift_pos_seq);
            l_pos_seq += (del_mask & shift_pos_seq);

            if (minz_pres.second){ // If at least one minimizer was present, search the kmer

                const const_UnitigColorMap<U> um = findUnitig(seqs.c_str(), pos_seq, seqs.length(), mhi);

                if (!um.isEmpty){

                    processUnitigMap(um, pos_seq);

                    pos_seq += um.len - 1;
                    mhi += pos_seq - mhi.getKmerPosition();
                }
            }

            ++pos_seq;
            ++mhi;

            while (pos_seq < seqs.length() - k_ + 1){

                shift_pos_seq = (pos_seq / k_) + (pos_seq % k_ > shift);
                l_pos_seq = pos_seq - (ins_mask & shift_pos_seq) + (del_mask & shift_pos_seq);

                it_min = *mhi;
                mhr = *it_min;

                // If minimizers of new kmer are different from minimizers of previous kmer
                // or if minimizers are the same but they were present, search them again
                if ((mhr.pos != minz_pres.first) || minz_pres.second){

                    if (mhr.pos != minz_pres.first){

                        minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                        minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                        for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                            mhr = *it_min;
                            minz = Minimizer(seqs.c_str() + mhr.pos).rep();
                            minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                        }
                    }

                    if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                        const const_UnitigColorMap<U> um = findUnitig(seqs.c_str(), pos_seq, seqs.length(), mhi);

                        if (!um.isEmpty){

                            processUnitigMap(um, pos_seq);

                            pos_seq += um.len - 1;
                            mhi += pos_seq - mhi.getKmerPosition();
                        }
                    }
                }

                ++pos_seq;
                ++mhi;
            }
        }
    };

    if (exact){

        for (size_t i = 0; i < seq.length() - k_ + 1; ++i) {

            const const_UnitigColorMap<U> um = findUnitig(seq.c_str(), i, seq.length());

            if (!um.isEmpty) { // Read maps to a Unitig

                if (um.strand){

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({i + j - um.dist, um.getKmerMapping(j)});
                }
                else {

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({i + um.dist + um.len - j - 1, um.getKmerMapping(j)});
                }

                i += um.len - 1;
            }
        }
    }

    if (substitution){

        for (size_t i = 0; i != k_; ++i){

            seqs = seq;
            worker_func(true, false, false, i);
        }
    }

    if (insertion){

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << seq[j];

            for (size_t j = i, cpt = 0; j < seq.length(); ++j, ++cpt) {

                if (cpt % (k_ - 1) == 0) ss << alpha[0];

                ss << seq[j];
            }

            seqs = ss.str();

            worker_func(false, true, false, i);
        }
    }

    if (deletion){

        const size_t sz = v_um.size();

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << seq[j];

            for (size_t j = i, cpt = 0; j < seq.length(); ++j, ++cpt) {

                if (cpt % (k_ + 1) != 0) ss << seq[j];
            }

            seqs = ss.str();

            worker_func(false, false, true, i);
        }
    }

    return v_um;
}*/

template<typename U>
bool ColoredCDBG<U>::search(const vector<string>& query_filenames, const string& out_filename_prefix,
                            const double ratio_kmers, const bool inexact_search, const size_t nb_threads) const {

    if (invalid){

        cerr << "ColoredCDBG::search(): Graph is invalid and cannot be searched" << endl;
        return false;
    }

    if (nb_threads > std::thread::hardware_concurrency()){

        cerr << "ColoredCDBG::search(): Number of threads cannot be greater than or equal to " << std::thread::hardware_concurrency() << "." << endl;
        return false;
    }

    if (nb_threads <= 0){

        cerr << "ColoredCDBG::search(): Number of threads cannot be less than or equal to 0." << endl;
        return false;
    }

    const string out_tmp = out_filename_prefix + ".tsv";

    FILE* fp_tmp = fopen(out_tmp.c_str(), "w");

    if (fp_tmp == NULL) {

        cerr << "ColoredCDBG::search(): Could not open file " << out_tmp << " for writing." << endl;
        return false;
    }
    else {

        fclose(fp_tmp);

        if (std::remove(out_tmp.c_str()) != 0) cerr << "ColoredCDBG::search(): Could not remove temporary file " << out_tmp << endl;
    }

    string s;

    size_t file_id = 0;

    const size_t k = this->getK();

    const size_t max_len_seq = 1024;
    const size_t thread_seq_buf_sz = 64 * max_len_seq;

    FileParser fp(query_filenames);

    ofstream outfile;
    ostream out(0);

    outfile.open(out_tmp.c_str());
    out.rdbuf(outfile.rdbuf());
    out.sync_with_stdio(false);

    const char query_pres[2] = {'\t', '1'};
    const char query_abs[2] = {'\t', '0'};

    const size_t l_query_res = 2;
    const size_t nb_colors = getNbColors();

    auto processCounts = [&](const vector<pair<size_t, const_UnitigColorMap<U>>>& v_um, Roaring* color_occ_r, uint32_t* color_occ_u){

        struct hash_pair {

            size_t operator() (const pair<size_t, pair<Kmer, size_t>>& p) const {

                return static_cast<size_t>(XXH64((const void *)&p, sizeof(pair<size_t, pair<Kmer, size_t>>), 0));
            }
        };

        unordered_set<pair<size_t, pair<Kmer, size_t>>, hash_pair> s_um;

        typename unordered_set<pair<size_t, pair<Kmer, size_t>>, hash_pair>::const_iterator it;

        for (const auto& p : v_um){

            s_um.insert({p.first, {p.second.strand ? p.second.getUnitigHead() : p.second.getUnitigTail().twin(), p.second.dist}});
        }

        for (const auto& p : v_um){

            const_UnitigColorMap<U> um = p.second;

            size_t pos_query = p.first;

            if (um.strand) {

                const Kmer head = um.getUnitigHead();

                size_t pos_unitig = um.dist;

                it = s_um.find({pos_query, {head, um.dist}});

                if (it != s_um.end()) {

                    s_um.erase(it);

                    while ((pos_unitig + k) < um.size){

                        ++pos_query;
                        ++pos_unitig;

                        it = s_um.find({pos_query, {head, pos_unitig}});

                        if (it == s_um.end()) break;
                        else {

                            ++(um.len);

                            s_um.erase(it);
                        }
                    }

                    const UnitigColors* uc = um.getData()->getUnitigColors(um);

                    UnitigColors::const_iterator it_uc = uc->begin(um);
                    UnitigColors::const_iterator it_uc_end = uc->end();

                    if (inexact_search){

                        for (; it_uc != it_uc_end; ++it_uc) color_occ_r[it_uc.getColorID()].add(it_uc.getKmerPosition() - um.dist + p.first);
                    }
                    else {

                        for (; it_uc != it_uc_end; ++it_uc) color_occ_u[it_uc.getColorID()] += 1;
                    }
                }
            }
            else {

                const Kmer head = um.getUnitigTail().twin();
                const size_t max_pos_um = um.dist + um.len - 1;

                it = s_um.find({pos_query, {head, um.dist}});

                if (it != s_um.end()) {

                    s_um.erase(it);

                    while (um.dist > 0){

                        ++pos_query;

                        it = s_um.find({pos_query, {head, um.dist - 1}});

                        if (it == s_um.end()) break;
                        else {

                            --(um.dist);
                            ++(um.len);

                            s_um.erase(it);
                        }
                    }

                    const UnitigColors* uc = um.getData()->getUnitigColors(um);

                    UnitigColors::const_iterator it_uc = uc->begin(um);
                    UnitigColors::const_iterator it_uc_end = uc->end();

                    if (inexact_search){

                        for (; it_uc != it_uc_end; ++it_uc) color_occ_r[it_uc.getColorID()].add(max_pos_um - it_uc.getKmerPosition() + p.first);
                    }
                    else {

                        for (; it_uc != it_uc_end; ++it_uc) color_occ_u[it_uc.getColorID()] += 1;
                    }
                }
            }
        }

        if (inexact_search){

            for (size_t i = 0; i < nb_colors; ++i) color_occ_u[i] = color_occ_r[i].cardinality();
        }
    };

    auto searchQuery = [&](const string& query, Roaring* color_occ_r, uint32_t* color_occ_u, const size_t nb_km_min){

        size_t nb_color_pres = 0;

        const vector<pair<size_t, const_UnitigColorMap<U>>> v_um_e = this->searchSequence(query, true, false, false, false, false);

        processCounts(v_um_e, color_occ_r, color_occ_u); // Extract k-mer occurrences for each color

        if (inexact_search){

            for (size_t j = 0, nb_color_pres = 0; j < nb_colors; ++j) nb_color_pres += (color_occ_u[j] >= nb_km_min);

            if (nb_color_pres == nb_colors) return;

            const vector<pair<size_t, const_UnitigColorMap<U>>> v_um_d = this->searchSequence(query, false, false, true, false, false);

            processCounts(v_um_d, color_occ_r, color_occ_u); // Extract k-mer occurrences for each color

            for (size_t j = 0, nb_color_pres = 0; j < nb_colors; ++j) nb_color_pres += (color_occ_u[j] >= nb_km_min);

            if (nb_color_pres == nb_colors) return;

            const vector<pair<size_t, const_UnitigColorMap<U>>> v_um_m = this->searchSequence(query, false, false, false, true, false);

            processCounts(v_um_m, color_occ_r, color_occ_u); // Extract k-mer occurrences for each color

            for (size_t j = 0, nb_color_pres = 0; j < nb_colors; ++j) nb_color_pres += (color_occ_u[j] >= nb_km_min);

            if (nb_color_pres == nb_colors) return;

            const vector<pair<size_t, const_UnitigColorMap<U>>> v_um_i = this->searchSequence(query, false, true, false, false, false);

            processCounts(v_um_i, color_occ_r, color_occ_u); // Extract k-mer occurrences for each color
        }
    };

    auto writeCounts = [&](const char* query_name, const uint32_t* color_occ, char* buffer_res, size_t& pos_buffer_out, const size_t nb_km_min){

        const size_t l_query_name = strlen(query_name);

        if (pos_buffer_out + l_query_name >= thread_seq_buf_sz){ // If next result cannot fit in the buffer

            out.write(buffer_res, pos_buffer_out); // Write result buffer
            pos_buffer_out = 0; // Reset position to 0;
        }

        // Copy new result to buffer
        std::memcpy(buffer_res + pos_buffer_out, query_name, l_query_name * sizeof(char));

        pos_buffer_out += l_query_name;

        for (size_t i = 0; i < nb_colors; ++i){

            if (pos_buffer_out + l_query_res >= thread_seq_buf_sz){ // If next result cannot fit in the buffer

                out.write(buffer_res, pos_buffer_out); // Write result buffer
                pos_buffer_out = 0; // Reset position to 0;
            }

            std::memcpy(buffer_res + pos_buffer_out, (color_occ[i] >= nb_km_min) ? query_pres : query_abs, l_query_res * sizeof(char));

            pos_buffer_out += l_query_res;
        }

        if (pos_buffer_out + 1 >= thread_seq_buf_sz){ // If next result cannot fit in the buffer

            out.write(buffer_res, pos_buffer_out); // Write result buffer
            pos_buffer_out = 0; // Reset position to 0;
        }

        buffer_res[pos_buffer_out] = '\n';

        ++pos_buffer_out;
    };

    // Write header to TSV file
    {
        const vector<string> color_names = getColorNames();

        out << "query_name";

        for (const auto& name : color_names) out << '\t' << name;

        out << '\n';
    }

    if (nb_threads == 1){

        char* buffer_res = new char[thread_seq_buf_sz];

        uint32_t* color_occ_u = new uint32_t[nb_colors]();

        Roaring* color_occ_r = inexact_search ? new Roaring[nb_colors] : nullptr;

        size_t pos_buffer_out = 0;

        while (fp.read(s, file_id)){

            const size_t nb_km_min = static_cast<double>(s.length() - k + 1) * ratio_kmers;

            searchQuery(s, color_occ_r, color_occ_u, nb_km_min);
            writeCounts(fp.getNameString(), color_occ_u, buffer_res, pos_buffer_out, nb_km_min); // Output presence/absence for each color

            std::memset(color_occ_u, 0, nb_colors * sizeof(uint32_t));

            if (inexact_search){

                for (size_t j = 0; j < nb_colors; ++j) color_occ_r[j] = Roaring(); // Reset k-mer occurences for each color
            }
        }

        // Flush unresult written to final output
        if (pos_buffer_out > 0) out.write(buffer_res, pos_buffer_out);

        delete[] buffer_res;
        delete[] color_occ_u;

        if (color_occ_r != nullptr) delete[] color_occ_r;
    }
    else {

        bool stop = false;

        vector<thread> workers; // need to keep track of threads so we can join them
        vector<vector<string>> buffers_seq(nb_threads);
        vector<vector<string>> buffers_name(nb_threads);

        char** buffer_res = new char*[nb_threads];
        uint32_t** color_occ_u = new uint32_t*[nb_threads];
        Roaring** color_occ_r = new Roaring*[nb_threads];

        mutex mutex_files_in, mutex_file_out;

        for (size_t t = 0; t < nb_threads; ++t){

            buffer_res[t] = new char[thread_seq_buf_sz];
            color_occ_u[t] = new uint32_t[nb_colors]();
            color_occ_r[t] = inexact_search ? new Roaring[nb_colors] : nullptr;

            workers.emplace_back(

                [&, t]{

                    while (true) {

                        {
                            if (stop) return;

                            size_t buffer_sz = 0;

                            unique_lock<mutex> lock(mutex_files_in);

                            stop = !fp.read(s, file_id);

                            while (!stop){

                                buffer_sz += s.length();

                                buffers_seq[t].push_back(move(s));
                                buffers_name[t].push_back(string(fp.getNameString()));

                                if (buffer_sz >= thread_seq_buf_sz) break;
                                else stop = !fp.read(s, file_id);
                            }
                        }

                        size_t pos_buffer_out = 0;

                        const size_t buffers_seq_sz = buffers_seq[t].size();

                        for (size_t i = 0; i < buffers_seq_sz; ++i){

                            const size_t nb_km_min = static_cast<double>(buffers_seq[t][i].length() - k + 1) * ratio_kmers;

                            searchQuery(buffers_seq[t][i], color_occ_r[t], color_occ_u[t], nb_km_min);

                            if ((pos_buffer_out + buffers_name[t][i].length() + nb_colors * l_query_res + 2) >= thread_seq_buf_sz){

                                unique_lock<mutex> lock(mutex_file_out);

                                writeCounts(buffers_name[t][i].c_str(), color_occ_u[t], buffer_res[t], pos_buffer_out, nb_km_min); // Output presence/absence for each color
                            }
                            else writeCounts(buffers_name[t][i].c_str(), color_occ_u[t], buffer_res[t], pos_buffer_out, nb_km_min);

                            std::memset(color_occ_u[t], 0, nb_colors * sizeof(uint32_t));

                            if (inexact_search){

                                for (size_t j = 0; j < nb_colors; ++j) color_occ_r[t][j] = Roaring(); // Reset k-mer occurences for each color
                            }
                        }

                        if (pos_buffer_out > 0){

                            unique_lock<mutex> lock(mutex_file_out);

                            out.write(buffer_res[t], pos_buffer_out);
                        }

                        // Clear buffers for next round
                        buffers_seq[t].clear();
                        buffers_name[t].clear();
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        for (size_t t = 0; t < nb_threads; ++t){

            delete[] buffer_res[t];
            delete[] color_occ_u[t];

            if (color_occ_r[t] != nullptr) delete[] color_occ_r[t];
        }

        delete[] buffer_res;
        delete[] color_occ_u;
        delete[] color_occ_r;
    }

    outfile.close();
    fp.close();

    return true;
}

template<typename U>
void ColoredCDBG<U>::checkColors(const vector<string>& filename_seq_in) const {

    cout << "ColoredCDBG::checkColors(): Start" << endl;

    size_t file_id = 0;

    string s;

    KmerHashTable<tiny_vector<size_t, 1>> km_h;

    FastqFile FQ(filename_seq_in);

    while (FQ.read_next(s, file_id) >= 0){

        for (KmerIterator it_km(s.c_str()), it_km_end; it_km != it_km_end; ++it_km) {

            pair<KmerHashTable<tiny_vector<size_t, 1>>::iterator, bool> it = km_h.insert(it_km->first.rep(), tiny_vector<size_t, 1>());

            tiny_vector<size_t, 1>& tv = *(it.first);

            const size_t id = file_id / 64;

            while (tv.size() < (id + 1)) tv.push_back(0);

            tv[id] |= (1ULL << (file_id % 64));
        }
    }

    FQ.close();

    cout << "ColoredCDBG::checkColors(): All k-mers in the hash table with their colors" << endl;

    for (typename KmerHashTable<tiny_vector<size_t, 1>>::const_iterator it_km = km_h.begin(), it_km_end = km_h.end(); it_km != it_km_end; ++it_km){

        const Kmer km = it_km.getKey();
        const const_UnitigColorMap<U> ucm = this->find(km);

        if (ucm.isEmpty){

            cerr << "ColoredCDBG::checkColors(): K-mer " << km.toString() << " is not found in the graph" << endl;
            exit(1);
        }

        const UnitigColors* cs = ucm.getData()->getUnitigColors(ucm);

        if (cs == nullptr){

            cerr << "ColoredCDBG::checkColors(): K-mer " << km.toString() << " has no color set associated" << endl;
            exit(1);
        }

        const tiny_vector<size_t, 1>& tv = *it_km;
        const size_t tv_nb_max_elem = tv.size() * 64;

        for (size_t i = 0; i < std::min(filename_seq_in.size(), tv_nb_max_elem); ++i){

            const bool color_pres_graph = cs->contains(ucm, i);
            const bool color_pres_hasht = ((tv[i/64] >> (i%64)) & 0x1) == 0x1;

            if (color_pres_graph != color_pres_hasht){

                cerr << "ColoredCDBG::checkColors(): Current color is " << i << ": " << filename_seq_in[i] << endl;
                cerr << "ColoredCDBG::checkColors(): K-mer " << km.toString() << " for color " << i << ": " << filename_seq_in[i] << endl;
                cerr << "ColoredCDBG::checkColors(): Size unitig: " << ucm.size << endl;
                cerr << "ColoredCDBG::checkColors(): Mapping position: " << ucm.dist << endl;
                cerr << "ColoredCDBG::checkColors(): Mapping strand: " << ucm.strand << endl;
                cerr << "ColoredCDBG::checkColors(): Present in graph: " << color_pres_graph << endl;
                cerr << "ColoredCDBG::checkColors(): Present in hash table: " << color_pres_hasht << endl;

                exit(1);
            }
        }
    }

    cout << "ColoredCDBG::checkColors(): Checked all colors of all k-mers: everything is fine" << endl;
    cout << "ColoredCDBG::checkColors(): Number of k-mers in the graph: " << km_h.size() << endl;
}

#endif
