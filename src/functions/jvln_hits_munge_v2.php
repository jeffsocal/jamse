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

function noralized_value($x, $n)
{
    return (($x * $x) / $n);
}

function jvln_hits_munge_v2(array $hits, array $peaks, float $pre_mass, float $pre_z, int $n_hits)
{
    $frg = new TandemHash('yb', '12', 'w');

    $hits = table_invert($hits);

    if (key_exists('factors', $hits)) {
        $factors = $hits['factors'];
        unset($hits['factors']);

        foreach ($factors as $i => $terms) {
            $terms = json_decode($terms, TRUE);
            // exit;
            $hits['bm15'][$i] = $terms['bm15'];
            $hits['bm25'][$i] = $terms['bm25a'];

            foreach ($terms['fields'][0] as $term => $val) {

                if (preg_match("/idf|wlccs|atc/", $term))
                    $hits[$term][$i] = $val;
            }
        }
    }

    $sum_scores = array_sum($hits['score']);
    $hits['score'] = array_map("noralized_value", $hits['score'], array_fill(0, sizeof($hits['score']), $sum_scores));
    $scores = array_map("log", $hits['score']);

    $s_mean = array_mean($scores);
    $s_sdev = max(0.001, array_stdev($scores));

    $lnorm = new Normal($s_mean, $s_sdev);

    $hits = table_head($hits, ceil($n_hits * 1.5));
    $hits = table_fillCol($hits, 'pre_mma', 0);
    $hits = table_fillCol($hits, 'psm_ovr_y', 0);
    $hits = table_fillCol($hits, 'psm_ovr_b', 0);
    $hits = table_fillCol($hits, 'psm_z', 0);
    $hits = table_fillCol($hits, 'psm_i', 0);
    $hits = table_fillCol($hits, 'psm_mma', 0);
    $hits = table_fillCol($hits, 'psm_score', 0);
    $hits = table_fillCol($hits, 'psm_fcorr', 0);
    $hits = table_fillCol($hits, 'psm_fcorr_y', 0);
    $hits = table_fillCol($hits, 'psm_fcorr_b', 0);
    $hits = table_fillCol($hits, 'hit_prob', 0);
    $hits = table_fillCol($hits, 'psm_prob', 0);
    $hits = table_fillCol($hits, 'alt_prob', 0);

    for ($i = 0; $i < table_length($hits); $i ++) {
        $score = $hits['score'][$i];
        $peptide = $hits['peptide'][$i];

        $z_table = [];
        for ($z = 0; $z <= 5; $z ++) {
            $this_mass = chargeMass($frg->getMolecularWeight($peptide), $z);
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

        $usr_z = '1';
        if ($mmz == 3)
            $usr_z = '12';

        $frg_y = new TandemHash('y', $usr_z, 'w');
        $frg_b = new TandemHash('b', $usr_z, 'w');
        $frg_a = new TandemHash('yb', $usr_z, 'w');

        $ovr_y = $frg_y->hashOverlap($peptide, $peaks);
        $ovr_b = $frg_b->hashOverlap($peptide, $peaks);
        $ovr_a = $frg_a->hashOverlap($peptide, $peaks);

        $frg_mean = array_mean($ovr_a['overlap']);
        $frg_stdv = array_stdev($ovr_a['overlap']);

        $hits['pre_mma'][$i] = round($mma, 4);
        $hits['psm_ovr_y'][$i] = count($ovr_y['overlap']);
        $hits['psm_ovr_b'][$i] = count($ovr_b['overlap']);
        $hits['psm_z'][$i] = $mmz;
        $hits['psm_i'][$i] = $mmi;
        $hits['psm_mma'][$i] = round($frg_mean, 4);

        /*
         * spectral properties
         */
        $n_mass = neutralMass($pre_mass, $mmz);

        $hits['psm_score'][$i] = sigfigs($score, 2);
        $hits['psm_fcorr'][$i] = sigfigs($ovr_a['overlap_score'], 3);
        $hits['psm_fcorr_y'][$i] = sigfigs($ovr_y['overlap_score'], 3);
        $hits['psm_fcorr_b'][$i] = sigfigs($ovr_b['overlap_score'], 3);
        $hits['hit_prob'][$i] = sigfigs(abs($lnorm->cdf(log($score))), 5);

        $hits['psm_prob'][$i] = '';
        $hits['alt_prob'][$i] = 0;
    }

    $hits = table_sort($hits, 'psm_score', 'desc', 'number');

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
                $hits['psm_prob'][$i] = 0;
                $hits['alt_prob'][$i] = sigfigs(bayesimple($prob_this, 1 - $prob_max), 3);
            }
            $drop_rows[] = $i;
            continue;
        }
        $hits['psm_prob'][$i] = sigfigs(bayesimple($prob_this, 1 - $prob_max), 3);
    }
    $alts = table_keeprows($hits, $drop_rows);
    $alts = table_droprows_lt($alts, 'alt_prob', 0.001);
    $hits = table_droprows($hits, $drop_rows);

    $hits = table_reindex($hits);
    $hits['psm_prob'][0] = 1;
    if (table_length($hits) > 1)
        $hits['psm_prob'][0] = sigfigs(bayesimple($hits['hit_prob'][0], 1 - $hits['hit_prob'][1]), 3);

    $hits['hit_rank'] = range(1, table_length($hits));
    $hits = table_dropcols($hits, [
        'hit_prob',
        'score'
    ]);
    $hits = table_head($hits, $n_hits);

    $out['hits'] = $hits;
    $out['alts'] = $alts;
    return $out;
}
?>