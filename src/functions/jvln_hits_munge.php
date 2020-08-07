<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

use function BenTools\CartesianProduct\cartesian_probroduct;
use Jvln\Engines\TandemHash;
use MathPHP\Probability\Distribution\Discrete\Poisson;
use MathPHP\Probability\Distribution\Continuous\Normal;
use MathPHP\Probability\Distribution\Continuous\LogNormal;

function jvln_hits_munge(array $hits, array $peaks, float $pre_mass, float $pre_z, int $n_hits)
{
    $frg = new TandemHash('yb', '12', 'w');

    $scores = [];
    foreach ($hits as $vals) {
        $scores[] = log($vals['score']);
    }
    $s_mean = array_mean($scores);
    $s_sdev = max(0.001, array_stdev($scores));
    $lnorm = new Normal($s_mean, $s_sdev);

    // echo "\n mean:".$s_mean."\n sdev:".$s_sdev."\n";
    /*
     * pair down the hit table save some time
     */
    $hits = array_slice($hits, 0, $n_hits * 2);
    foreach ($hits as $row => $vals) {

        unset($hits[$row]['score']);

        $score = $vals['score'];
        $this_probeptide = $hits[$row]['peptide'];
        if(is_null($this_probeptide)){
            unset($hits[$row]);
            continue;
        }
        
        $this_overlap = $frg->hashOverlap($this_probeptide, $peaks);

        $frg_mean = array_mean($this_overlap['overlap']);
        $frg_stdv = array_stdev($this_overlap['overlap']);

        // foreach ($this_overlap as $i => $val) {
        // if (abs($val - $frg_mean) > $frg_stdv * 3)
        // unset($this_overlap[$i]);
        // }

        $z_table = [];
        for ($z = 0; $z <= 5; $z ++) {
            $this_mass = chargeMass($frg->getMolecularWeight($this_probeptide), $z);
            $z_table['z'][] = $z;
            $z_table['abs'][] = abs($pre_mass - $this_mass);
            $z_table['mma'][] = $pre_mass - $this_mass;
        }

        $z_table = table_sort($z_table, 'abs', 'ascn', 'number');

        $mma = $z_table['mma'][0];
        $mmz = $z_table['z'][0];
        $mmi = round(abs($mma * $mmz));

        if ($mmi >= 1)
            $mma = $pre_mass - chargeMass($frg->getMolecularWeight($this_probeptide) - $mmi, $mmz);

        $hits[$row]['pre_mma'] = round($mma, 4);
        $hits[$row]['psm_o'] = count($this_overlap['overlap']);
        $hits[$row]['psm_z'] = $mmz;
        $hits[$row]['psm_i'] = $mmi;
        // $hits[$row]['psm_r'] = round(count($this_overlap) / count($peaks), 3);
        $hits[$row]['psm_mma'] = round($frg_mean, 4);

        /*
         * spectral properties
         */
        $n_mass = neutralMass($pre_mass, $mmz);

        // /*
        // * information content probability
        // * log-normal distribution
        // */
        // $p_mean = 1.1420034;
        // $p_sd = 0.2178803;
        // $psm_lnorm = new LogNormal($p_mean, $p_sd);
        // $pvalue_info = $psm_lnorm->cdf($score / count($this_overlap['overlap']));
        // if ($pvalue_info . "" == "NAN")
        // $pvalue_info = 0.00;

        // /*
        // * peak overlap probability
        // * log-normal distribution
        // * fit from poly^2 :: overlap ~ pre_nmass from large dataset
        // */
        // $p_mean = - 0.5577285 * pow(log($n_mass), 2) + 8.6588187 * log($n_mass) - 31.8397488;
        // $p_sd = 0.06632917 * pow(log($n_mass), 2) - 1.17749888 * log($n_mass) + 5.43862010;
        // $psm_lnorm = new LogNormal($p_mean, $p_sd);
        // $pvalue_over = $psm_lnorm->cdf(count($this_overlap['overlap']));
        // if ($pvalue_over . "" == "NAN")
        // $pvalue_over = 0.00;

        /*
         * information recall probability
         * log-normal distribution
         * fit from poly^2 :: overlap ~ pre_nmass from large dataset
         */
        $p_mean = pow(log($n_mass), 2) * - 0.5920383 + log($n_mass) * 8.7718279 - 29.7646591;
        $p_sd = pow(log($n_mass), 2) * - 0.06369475 + log($n_mass) * 0.86396084 - 2.46828955;
        $psm_lnorm = new LogNormal($p_mean, $p_sd);
        $pvalue_recall = $psm_lnorm->cdf($score / 10);
        if ($pvalue_recall . "" == "NAN")
            $pvalue_recall = 0.00;

        $hits[$row]['psm_score'] = sigfigs($score / 10, 2);
        $hits[$row]['psm_fcorr'] = sigfigs($this_overlap['overlap_score'], 2);
        $hits[$row]['psm_prob'] = sigfigs(abs($pvalue_recall), 5);
        $dbs_prob = sigfigs(abs($lnorm->cdf(log($score))), 4);
        $hits[$row]['hit_prob'] = sigfigs(($hits[$row]['psm_prob'] * $dbs_prob), 5);

        $hits[$row]['rank_prob'] = '';
        $hits[$row]['alt_prob'] = 0;
    }

    $hits = table_invert($hits);
    $hits = table_sort($hits, 'psm_prob', 'desc', 'number');

    $peptide_r1 = getRootSequence($hits['peptide'][0]);

    if (strlen($peptide_r1) > 250) {
        echo "\n--WARN--\tpeptide.length: " . $peptide_r1 . "\n";
        $peptide_r1 = substr($peptide_r1, 1, 250);
    }

    $prob_max = $hits['hit_prob'][0];
    $alt_id = '';
    $alt_pb = 0;
    $drop_rows = [];
    for ($i = 1; $i < table_length($hits); $i ++) {
        $peptide_rn = getRootSequence($hits['peptide'][$i]);
        
        if (strlen($peptide_rn) > 250) {
            echo "\n--WARN--\tpeptide.length: " . $peptide_rn . "\n";
            $peptide_rn = substr($peptide_rn, 1, 250);
        }
        
        $prob_this = $hits['hit_prob'][$i];
        if (levenshtein($peptide_r1, $peptide_rn) <= 2) {
            if ($alt_id == '') {
                $hits['rank_prob'][$i] = 0;
                $hits['alt_prob'][$i] = sigfigs(bayesimple($prob_this, 1 - $prob_max), 3);
            }
            $drop_rows[] = $i;
            continue;
        }
        $hits['rank_prob'][$i] = sigfigs(bayesimple($prob_this, 1 - $prob_max), 3);
    }
    $alts = table_keeprows($hits, $drop_rows);
    $alts = table_droprows_lt($alts, 'alt_prob', 0.001);
    $hits = table_droprows($hits, $drop_rows);

    $hits = table_reindex($hits);
    $hits['rank_prob'][0] = 1;
    if (table_length($hits) > 1)
        $hits['rank_prob'][0] = sigfigs(bayesimple($hits['hit_prob'][0], 1 - $hits['hit_prob'][1]), 3);

    $hits['hit_rank'] = range(1, table_length($hits));
    $hits = table_head($hits, $n_hits);

    $out['hits'] = $hits;
    $out['alts'] = $alts;
    return $out;
}
?>