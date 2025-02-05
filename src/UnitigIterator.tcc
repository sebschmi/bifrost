#ifndef BIFROST_UNITIG_ITERATOR_TCC
#define BIFROST_UNITIG_ITERATOR_TCC

#include "CompactedDBG.hpp"

template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const>::unitigIterator() :  i(0), v_unitigs_sz(0), v_kmers_sz(0), h_kmers_ccov_sz(0), sz(0), invalid(true), cdbg(nullptr) {}

template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const>::unitigIterator(CompactedDBG_ptr_t cdbg_) :
                i(0), v_unitigs_sz(0), v_kmers_sz(0), h_kmers_ccov_sz(0), sz(0), invalid(true), cdbg(cdbg_),
                it_h_kmers_ccov((cdbg_ == nullptr) || cdbg_->invalid ? typename KmerHashTable<CompressedCoverage_t<U>>::const_iterator() : cdbg_->h_kmers_ccov.begin()){

    if ((cdbg != nullptr) && !cdbg->invalid && (cdbg->size() != 0)){

        invalid = false;

        v_unitigs_sz = cdbg->v_unitigs.size();
        v_kmers_sz = cdbg->km_unitigs.size();
        h_kmers_ccov_sz = cdbg->h_kmers_ccov.size();

        sz = v_unitigs_sz + v_kmers_sz + h_kmers_ccov_sz;
    }
}

template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const>::unitigIterator(const unitigIterator& o) :   i(o.i), v_unitigs_sz(o.v_unitigs_sz), v_kmers_sz(o.v_kmers_sz),
                                                                            it_h_kmers_ccov(o.it_h_kmers_ccov), h_kmers_ccov_sz(o.h_kmers_ccov_sz),
                                                                            sz(o.sz), invalid(o.invalid), um(o.um), cdbg(o.cdbg) {}

template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const>& unitigIterator<U, G, is_const>::operator++() {

    if (invalid) return *this;

    if ((cdbg == nullptr) || cdbg->invalid || (i >= sz)){

        invalid = true;
        return *this;
    }

    if (i < v_unitigs_sz){

        um = UnitigMap<U, G, is_const>(i, 0, cdbg->v_unitigs[i]->getSeq().size() - cdbg->getK() + 1,
                                       cdbg->v_unitigs[i]->getSeq().size(), false, false, true, cdbg);
    }
    else if (i < (v_unitigs_sz + v_kmers_sz)){

        um = UnitigMap<U, G, is_const>(i - v_unitigs_sz, 0, 1, cdbg->getK(), true, false, true, cdbg);
    }
    else {

        um = UnitigMap<U, G, is_const>(it_h_kmers_ccov.getHash(), 0, 1, cdbg->getK(), false, true, true, cdbg);

        ++it_h_kmers_ccov;
    }

    ++i;

    return *this;
}

template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const> unitigIterator<U, G, is_const>::operator++(int) {

    unitigIterator<U, G, is_const> tmp(*this);
    operator++();

    return tmp;
}

template<typename U, typename G, bool is_const>
bool unitigIterator<U, G, is_const>::operator==(const unitigIterator& o) const {

    if (invalid || o.invalid) return invalid && o.invalid;
    return  (i == o.i) && (v_unitigs_sz == o.v_unitigs_sz) && (v_kmers_sz == o.v_kmers_sz) &&
            (h_kmers_ccov_sz == o.h_kmers_ccov_sz) && (sz == o.sz) && (it_h_kmers_ccov == o.it_h_kmers_ccov) &&
            (cdbg == o.cdbg) && (um == o.um);
}

template<typename U, typename G, bool is_const>
bool unitigIterator<U, G, is_const>::operator!=(const unitigIterator& o) const { return !operator==(o); }

template<typename U, typename G, bool is_const>
const UnitigMap<U, G, is_const>& unitigIterator<U, G, is_const>::operator*() const { return um; }

template<typename U, typename G, bool is_const>
const UnitigMap<U, G, is_const>* unitigIterator<U, G, is_const>::operator->() const { return &um; }

template<typename U, typename G, bool is_const>
void unitigIterator<U, G, is_const>::setIndex(size_t i) {
    this->i = i;

    if (invalid) return;

    if ((cdbg == nullptr) || cdbg->invalid || (i >= sz)){

        invalid = true;
        return;
    }

    if (i < v_unitigs_sz){

        um = UnitigMap<U, G, is_const>(i, 0, cdbg->v_unitigs[i]->getSeq().size() - cdbg->getK() + 1,
                                       cdbg->v_unitigs[i]->getSeq().size(), false, false, true, cdbg);
    }
    else if (i < (v_unitigs_sz + v_kmers_sz)){

        um = UnitigMap<U, G, is_const>(i - v_unitigs_sz, 0, 1, cdbg->getK(), true, false, true, cdbg);
    }
    else {

        it_h_kmers_ccov = cdbg->h_kmers_ccov.begin() + (i - (v_unitigs_sz + v_kmers_sz));
        um = UnitigMap<U, G, is_const>(it_h_kmers_ccov.getHash(), 0, 1, cdbg->getK(), false, true, true, cdbg);
    }
}

template<typename U, typename G, bool is_const>
void unitigIterator<U, G, is_const>::setIndex(size_t i, const vector<size_t>& h_kmer_ccov_orders) {
    this->i = i;

    if (invalid) return;

    if ((cdbg == nullptr) || cdbg->invalid || (i >= sz)){

        invalid = true;
        return;
    }

    if (i < v_unitigs_sz){

        um = UnitigMap<U, G, is_const>(i, 0, cdbg->v_unitigs[i]->getSeq().size() - cdbg->getK() + 1,
                                       cdbg->v_unitigs[i]->getSeq().size(), false, false, true, cdbg);
    }
    else if (i < (v_unitigs_sz + v_kmers_sz)){

        um = UnitigMap<U, G, is_const>(i - v_unitigs_sz, 0, 1, cdbg->getK(), true, false, true, cdbg);
    }
    else {

        it_h_kmers_ccov = cdbg->h_kmers_ccov.begin();
        auto local_i = (i - (v_unitigs_sz + v_kmers_sz));
        it_h_kmers_ccov.h = h_kmer_ccov_orders[local_i];
        um = UnitigMap<U, G, is_const>(it_h_kmers_ccov.getHash(), 0, 1, cdbg->getK(), false, true, true, cdbg);
    }
}

template<typename U, typename G, bool is_const>
std::ostream& operator<<(std::ostream& stream, const unitigIterator<U, G, is_const>& it) {
    stream << "UnitigIterator(i = " << it.i;
    stream << ", v_unitigs_sz = " << it.v_unitigs_sz;
    stream << ", v_kmers_sz = " << it.v_kmers_sz;
    stream << ", h_kmers_ccov_sz = "<< it.h_kmers_ccov_sz;
    stream << ", sz = " << it.sz;
    stream << ", invalid = " << it.invalid;
    stream << ", it_h_kmers_ccov = ";
    stream.flush();
    append_kmer_hash_table_iterator_to_stream<CompressedCoverage_t<U>, true>(stream, it.it_h_kmers_ccov);
    stream.flush();
    stream << ", um = ";
    append_unitig_map_to_stream_robust(stream, it.um);
    stream.flush();
    stream << ", cdbg = " << (void*) it.cdbg << ")";
    stream.flush();
    return stream;
}

#endif
