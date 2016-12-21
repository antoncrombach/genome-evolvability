//
// Mutating parameters can be done in different ways.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "mutate_rates.hh"
#include "chromosome.hh"


void
fluke::FixedMutateRates::mutate( Chromosome &c ) {
    if( uniform() < 1.0 - pow( 1.0 - c.mut_rate_, c.NR_RATES ) ) {
        uint aux = rand_range( c.NR_RATES );
        switch( aux ) {
        case 0: c.cp_gene_rate_ += c.mut_step_;
            break;
        case 1: c.rm_gene_rate_ += c.mut_step_;
            break;
        case 2: c.dsb_recombination_ += c.dsb_step_;
            break;
        case 3: c.cp_gene_rate_ -= c.mut_step_;
            if( c.cp_gene_rate_ < 0.0 ) c.cp_gene_rate_ = 0.0;
            break;
        case 4: c.rm_gene_rate_ -= c.mut_step_;
            if( c.rm_gene_rate_ < 0.0 ) c.rm_gene_rate_ = 0.0;
            break;
        case 5: c.dsb_recombination_ -= c.dsb_step_;
            if( c.dsb_recombination_ < 0.0 ) c.dsb_recombination_ = 0.0;
            break;
        // and some extra cases for retroposon activity
        case 6: c.cp_tp_rate_ += c.retro_step_;
            break;
        case 7: c.rm_tp_rate_ += c.retro_step_;
            break;
        case 8: c.rm_ltr_rate_ += c.retro_step_;
            break;
        case 9: c.cp_tp_rate_ -= c.retro_step_;
            if( c.cp_tp_rate_ < 0.0 ) c.cp_tp_rate_ = 0.0;
            break;
        case 10: c.rm_tp_rate_ -= c.retro_step_;
            if( c.rm_tp_rate_ < 0.0 ) c.rm_tp_rate_ = 0.0;
            break;
        case 11: c.rm_ltr_rate_ -= c.retro_step_;
            if( c.rm_ltr_rate_ < 0.0 ) c.rm_ltr_rate_ = 0.0;
            break;
        }
    }
}

void
fluke::LinearMutateRates::mutate( Chromosome &c ) {
    // NOTE: if rate hits zero, it stays zero!
    // Maybe introduce minimal mutation rate?
    if( uniform() < 1.0 - pow( 1.0 - c.mut_rate_, c.NR_RATES ) ) {
        uint aux = rand_range( c.NR_RATES );
        switch( aux ) {
        case 0: c.cp_gene_rate_ += c.mut_step_ * c.cp_gene_rate_;
            break;
        case 1: c.rm_gene_rate_ += c.mut_step_ * c.rm_gene_rate_;
            break;
        case 2: c.dsb_recombination_ += c.dsb_step_;
            break;
        case 3: c.cp_gene_rate_ -= c.mut_step_ * c.cp_gene_rate_;
            if( c.cp_gene_rate_ < 0.0 ) c.cp_gene_rate_ = 0.0;
            break;
        case 4: c.rm_gene_rate_ -= c.mut_step_ * c.rm_gene_rate_;
            if( c.rm_gene_rate_ < 0.0 ) c.rm_gene_rate_ = 0.0;
            break;
        case 5: c.dsb_recombination_ -= c.dsb_step_;
            if( c.dsb_recombination_ < 0.0 ) c.dsb_recombination_ = 0.0;
            break;
        // and some extra cases for retroposon activity
        case 6: c.cp_tp_rate_ += c.retro_step_;
            break;
        case 7: c.rm_tp_rate_ += c.retro_step_;
            break;
        case 8: c.rm_ltr_rate_ += c.retro_step_;
            break;
        case 9: c.cp_tp_rate_ -= c.retro_step_;
            if( c.cp_tp_rate_ < 0.0 ) c.cp_tp_rate_ = 0.0;
            break;
        case 10: c.rm_tp_rate_ -= c.retro_step_;
            if( c.rm_tp_rate_ < 0.0 ) c.rm_tp_rate_ = 0.0;
            break;
        case 11: c.rm_ltr_rate_ -= c.retro_step_;
            if( c.rm_ltr_rate_ < 0.0 ) c.rm_ltr_rate_ = 0.0;
            break;
        }
    }
}

void
fluke::UniformMutateRates::mutate( Chromosome &c ) {
    // FIX How to do this????
    if( uniform() < 1.0 - pow( 1.0 - c.mut_rate_, c.NR_RATES ) ) {
        uint aux = rand_range( c.NR_RATES );
        double bux = rand_range( high_ - low_ ) + low_;
        // mutate
        switch( aux ) {
        case 0: c.cp_gene_rate_ = bux;
            break;
        case 1: c.rm_gene_rate_ = bux;
            break;
        case 2: c.dsb_recombination_ += c.dsb_step_;
            break;
        case 3: c.cp_gene_rate_ = bux;
            break;
        case 4: c.rm_gene_rate_ = bux;
            break;
        case 5: c.dsb_recombination_ -= c.dsb_step_;
            if( c.dsb_recombination_ < 0.0 ) c.dsb_recombination_ = 0.0;
            break;
        // and some extra cases for retroposon activity
        case 6: c.cp_tp_rate_ += c.retro_step_;
            break;
        case 7: c.rm_tp_rate_ += c.retro_step_;
            break;
        case 8: c.rm_ltr_rate_ += c.retro_step_;
            break;
        case 9: c.cp_tp_rate_ -= c.retro_step_;
            if( c.cp_tp_rate_ < 0.0 ) c.cp_tp_rate_ = 0.0;
            break;
        case 10: c.rm_tp_rate_ -= c.retro_step_;
            if( c.rm_tp_rate_ < 0.0 ) c.rm_tp_rate_ = 0.0;
            break;
        case 11: c.rm_ltr_rate_ -= c.retro_step_;
            if( c.rm_ltr_rate_ < 0.0 ) c.rm_ltr_rate_ = 0.0;
            break;
        }
    }
}
