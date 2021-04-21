#ifndef BIFROST_HELSITIGS_TCC
#define BIFROST_HELSITIGS_TCC

template<typename U, typename G>
bool CompactedDBG<U, G>::convert_tigs(CompactedDBG<U, G> dbg, const Tigs tigs, const size_t nb_threads) {
    cout << "convert_tigs(dbg, tigs, nb_threads = " << nb_threads << ")" << endl;

    cout << "given unitigs:" << endl;
    vector<string> unitigs;
    for (const auto unitig : dbg) {
        const string seq(unitig.referenceUnitigToString());
        cout << seq << endl;
        unitigs.push_back(seq);
    }

    hmap_min_unitigs = move(dbg.hmap_min_unitigs);
    for (const auto& unitig : unitigs) {
        addUnitig(unitig, (unitig.length() == k_) ? km_unitigs.size() : v_unitigs.size());
    }
}

#endif