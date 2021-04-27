#ifndef BIFROST_COLORED_HELSITIGS_TCC
#define BIFROST_COLORED_HELSITIGS_TCC

//#include<set>

template<>
void CompactedDBG<DataAccessor<void>, DataStorage<void>>::colorUnitig(const CompactedDBG<DataAccessor<void>,DataStorage<void>>* dbg,
                                                                      const UnitigMap<DataAccessor<void>, DataStorage<void>>& um,
                                                                      const ptrdiff_t* tigs_edge_out_offset,
                                                                      const ptrdiff_t* tigs_edge_out_limit,
                                                                      const size_t* tigs_insert_out_offset,
                                                                      const size_t* tigs_insert_out_limit,
                                                                      const size_t nb_threads,
                                                                      const bool initialise_data,
                                                                      const size_t total_tigs) {
    //cout << "\n === CompactedDBG::colorUnitig(): start for " << um.mappedSequenceToString() << endl;

    auto* self = dynamic_cast<ColoredCDBG<void>*>(this);
    if (self == 0) {
        cerr << "Could not cast this down to colored dbg" << endl;
        exit(1);
    }
    const auto* colored_dbg = dynamic_cast<const ColoredCDBG<void>*>(dbg);
    if (colored_dbg == 0) {
        cerr << "Could not cast dbg down to colored dbg" << endl;
        exit(1);
    }
    if (initialise_data) {
        DataStorage<void> new_ds(31, total_tigs, colored_dbg->getData()->color_names);
        *(self->getData()) = move(new_ds);
    }

    auto tig_edge = tigs_edge_out_offset;
    auto tig_insert = tigs_insert_out_offset;
    auto unitig_iterator = colored_dbg->begin();
    DataStorage<void>* ds = self->getData();
    //cout << "UnitigColors is " << (void*) ds->getUnitigColors(um) << endl;
    size_t offset = 0;

    *(um.getData()) = ds->insert(um).first;

    //set<pair<pair<size_t, bool>, size_t>> inserted_colors;

    for (; tig_edge < tigs_edge_out_limit && tig_insert < tigs_insert_out_limit; tig_edge++, tig_insert++) {
        auto edge = *tig_edge;
        auto insert = *tig_insert;

        if (insert == 0) {
            unitig_iterator.setIndex(abs(edge));
            const auto unitig_mapping = *unitig_iterator;
            //cout << "\nCopying colors from " << unitig_mapping.mappedSequenceToString() << endl;

            // process existing unitig
            // TODO process runs of colors at once
            const UnitigColors* unitig_colors = unitig_mapping.getData()->getUnitigColors(unitig_mapping);
            for (auto position_color = unitig_colors->begin(unitig_mapping); position_color != unitig_colors->end(); position_color++) {
                const auto position = (*position_color).first;
                const auto color = (*position_color).second;

                auto tmp_um = um;
                if (edge >= 0) {
                    tmp_um.dist = offset + position;
                    tmp_um.len = 1;
                } else {
                    tmp_um.dist = offset + position;
                    tmp_um.len = 1;
                    tmp_um.strand = !tmp_um.strand;
                }
                //cout << "Writing color " << color << " to strand " << tmp_um.strand << " of " << tmp_um.getIndex() << (tmp_um.isAbundant ? "a" : "") << (tmp_um.isShort ? "s" : "") << "[" << tmp_um.dist << "/" << tmp_um.size << "] = " << tmp_um.getMappedHead().toString() << endl;
                //ds->getUnitigColors(tmp_um)->add(tmp_um, color);
                //cout << "Before unitig colors: " << (void*) tmp_um.getData()->getUnitigColors(tmp_um) << ", Unitig data: " << (void*) tmp_um.getData() << ", UnitigMap cdbg: " << (void*) tmp_um.cdbg << ", pos_unitig: " << tmp_um.pos_unitig << ", flag = " << tmp_um.getData()->getUnitigColors(tmp_um)->flag() << endl;
                tmp_um.getData()->getUnitigColors(tmp_um)->add(tmp_um, color);
                //cout << " After unitig colors: " << (void*) tmp_um.getData()->getUnitigColors(tmp_um) << ", Unitig data: " << (void*) tmp_um.getData() << ", UnitigMap cdbg: " << (void*) tmp_um.cdbg << ", pos_unitig: " << tmp_um.pos_unitig << ", flag = " << tmp_um.getData()->getUnitigColors(tmp_um)->flag() << endl;
                //inserted_colors.insert({{tmp_um.dist, true}, color});
                //inserted_colors.insert({{tmp_um.dist, false}, color});

                /*if (!ds->getUnitigColors(tmp_um)->contains(tmp_um, color)) {
                    cout << "Writing color failed" << endl;
                    exit(1);
                }*/

                /*for (size_t dist = 0; dist < um.size - k_; dist++) {
                    for (size_t color = 0; color < colored_dbg->getNbColors(); color++) {
                        tmp_um.dist = dist;
                        tmp_um.strand = true;

                        bool stored = ds->getUnitigColors(tmp_um)->contains(tmp_um, color);
                        bool inserted = (inserted_colors.find({{dist, true}, color}) != inserted_colors.end());
                        if (stored != inserted) {
                            cout << "Written colors mismatch with stored colors" << endl;
                            cout << "dist = " << dist << "; strand = " << true << "; color = " << color << endl;
                            cout << "stored = " << stored << "; inserted = " << inserted << endl;
                            exit(1);
                        }

                        tmp_um.dist = dist;
                        tmp_um.strand = false;

                        stored = ds->getUnitigColors(tmp_um)->contains(tmp_um, color);
                        inserted = (inserted_colors.find({{dist, false}, color}) != inserted_colors.end());
                        if (stored != inserted) {
                            cout << "Written colors mismatch with stored colors" << endl;
                            cout << "dist = " << dist << "; strand = " << false << "; color = " << color << endl;
                            cout << "stored = " << stored << "; inserted = " << inserted << endl;
                            exit(1);
                        }
                    }
                }*/

                /*for (const auto& inserted_color : inserted_colors) {
                    tmp_um.dist = inserted_color.first.first;
                    tmp_um.strand = inserted_color.first.second;
                    if (!ds->getUnitigColors(tmp_um)->contains(tmp_um, inserted_color.second)) {
                        cerr << "Formerly written color disappeared or writing color failed" << endl;
                        exit(1);
                    }
                }*/
            }

            offset += unitig_mapping.len;
        } else {
            // process dummy edge
            //cout << "Processing dummy edge" << endl;

            for (size_t position = 0; position < insert; position++) {
                const Kmer km = um.getMappedKmer(offset + position);
                const const_UnitigColorMap<void> source_unitig_mapping = colored_dbg->find(km);

                for (size_t color = 0; color < colored_dbg->getNbColors(); color++) {
                    if (source_unitig_mapping.getData()->getUnitigColors(source_unitig_mapping)->contains(source_unitig_mapping, color)) {
                        auto tmp_um = um;
                        if (edge >= 0) {
                            tmp_um.dist = offset + position;
                            tmp_um.len = 1;
                        } else {
                            tmp_um.dist = offset + position;
                            tmp_um.len = 1;
                            tmp_um.strand = !tmp_um.strand;
                        }

                        //cout << "Writing color " << color << " to strand " << tmp_um.strand << " of " << tmp_um.getIndex() << (tmp_um.isAbundant ? "a" : "") << (tmp_um.isShort ? "s" : "") << "[" << tmp_um.dist << "/" << tmp_um.size << "] = " << tmp_um.getMappedHead().toString() << endl;
                        tmp_um.getData()->getUnitigColors(tmp_um)->add(tmp_um, color);
                        //inserted_colors.insert({{tmp_um.dist, true}, color});
                        //inserted_colors.insert({{tmp_um.dist, false}, color});
                    }
                }
            }

            offset += insert;
        }
    }

    um.getData()->getUnitigColors(um)->optimizeFullColors(um);
}

#endif