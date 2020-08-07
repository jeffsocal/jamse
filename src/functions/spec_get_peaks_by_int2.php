<?php
use MathPHP\Statistics\Descriptive;
use MathPHP\Probability\Distribution\Continuous\Normal;

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

function spec_prob(float $spec_fcorr, int $spec_mz, int $spec_z)
{
    $spec_nmass = neutralMass($spec_mz, $spec_z);

    switch ($spec_z) {
        case 1:
            $mean = [
                - 3.631228e-01,
                5.500472e-03,
                - 2.558244e-06
            ];
            $sdev = [
                2.045233e-01,
                1.356978e-03,
                - 5.631944e-07
            ];
            break;
        case 2:
            $mean = [
                - 1.988416e+00,
                6.444450e-03,
                - 1.521552e-06
            ];
            $sdev = [
                4.618133e-01,
                7.647404e-04,
                - 1.887728e-07
            ];
            break;
        default:
            $mean = [
                - 1.491842e-01,
                2.511513e-03,
                - 4.021881e-07
            ];
            $sdev = [
                6.505080e-01,
                3.110282e-04,
                - 7.492960e-08
            ];
    }

    $p_mean = max(0.01, $mean[2] * pow(($spec_nmass), 2) + $mean[1] * ($spec_nmass) + $mean[0]);
    $p_sdev = max(0.01, $sdev[2] * pow(($spec_nmass), 2) + $sdev[1] * ($spec_nmass) + $sdev[0]);
    $prob = new Normal($p_mean, $p_sdev);
    return $prob->cdf($spec_fcorr);
}

function spec_get_peaks_by_int2($table, $pre_mz, $pre_z, $npeaks = 40, $difpeaks = 17.5)
{
    $amass = averagineMass();
    $mth = new Descriptive();

    $spec = $table;

    $int_sum = array_sum($spec['int']) * 1;
    $int_max = array_max($spec['int']) * 1;
    $int_mean = array_mean($spec['int']) * 1;
    $int_median = array_median($spec['int']) * 1;
    $int_sdev = array_stdev($spec['int']) * 1;
    $int_adev = array_absdev($spec['int']) * 1;
    $int_perc_10 = $mth->percentile($spec['int'], 10) * 1;
    $int_iqr = $mth->iqr($spec['int']) * 1;

    /*
     * normalized spectral abundance metrics
     * multiply by 1 to remove internal commas
     */
    $out['int_rsd'] = sigfigs($int_sdev / $int_mean, 3);
    $out['int_rad'] = sigfigs($int_adev / $int_median, 3);
    $out['int_iqr'] = sigfigs($int_iqr / $int_max, 3);
    $out['int_ntp'] = sigfigs($int_perc_10 / $int_max, 3);

    /*
     * remove an over abundance of peaks in the spectra
     */
    $spec = table_sort($spec, 'int', 'desc', 'number');
    $spec = table_head($spec, round($npeaks * 15));

    $spec['peak'] = range(0, table_length($spec) - 1, 1);
    $spec = $ref = array_rowtocol($spec);

    $pre_int = [
        0
    ];
    /*
     * cut out precursor peaks
     * deleting peaks that are within 2 Da of precuror
     */
    foreach ($spec as $i => $v) {
        if (abs($v['mz'] - $pre_mz) <= 2) {
            if (abs($v['mz'] - $pre_mz) <= 0.1)
                $pre_int[] = $v['int'];
            unset($spec[$i]);
        }
    }
    /*
     * record a value indicating partial fragmentation
     * multiply by 1 to remove internal commas
     */
    $out['int_prd'] = sigfigs(array_max($pre_int) / $int_median, 3);

    $spec = array_values($spec);

    /*
     * iterate through array removing peaks around
     * the tallest peaks by intensity
     */
    $toss = [];
    $keep = [];
    while (count($keep) < $npeaks) {

        if (count($spec) == 0)
            break;

        $mz = $spec[0]['mz'];
        $in = $spec[0]['int'];
        $pk = $spec[0]['peak'];
        $keep[] = $pk;

        unset($spec[0]);

        foreach ($spec as $i => $v) {
            $diff = $mz - $v['mz'];
            if (abs($diff) <= $difpeaks) {
                unset($spec[$i]);
                if (abs($diff - 1) <= 1.05 and $v['int'] > $in / 2) {

                    $toss[] = $pk;
                    $keep[] = $v['peak'];
                }
            }
        }
        $spec = array_values($spec);
    }

    $keep = array_diff($keep, $toss);

    $keep = array_intersect_key($ref, array_flip($keep));

    if (count($keep) == 0)
        $keep = $ref;

    $t_bins = averagineBins($pre_mz, 2);

    $bins = array_fill(0, $t_bins, 0);
    foreach ($keep as $i => $vals) {
        if ($vals['int'] < $int_median)
            continue;

        $bin = min($t_bins - 1, ceil($vals['mz'] / $amass) * 1);
        $bins[$bin] = min(2, $bins[$bin] + 1);
    }

    for ($i = 0; $i < $t_bins; $i ++) {
        $f_bins[$i] = $bins[$i] * ($i / ($t_bins - 1));
    }

    $keep = table_invert($keep);

    $out['spec_fcorr'] = sigfigs(array_sum($f_bins), 4);
    $out['spec_prob'] = sigfigs(spec_prob($out['spec_fcorr'], $pre_mz, $pre_z), 4);
    $out['t_peaks'] = count($table['mz']);
    $out['n_peaks'] = count($keep['mz']);
    $out['peaks'] = array_values($keep['mz']);

    return $out;
}

?>