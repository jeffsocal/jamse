<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace Jvln\Engines;

use ProxySci\Omics\Proteomics\Fragments;

class TandemHash extends Fragments
{

    private $precursor_mass_bin;

    private $frag_mass_max;

    private $frag_mass_min;

    private $frag_mass_cut;

    private $m_int;

    private $m_slp;

    function __construct($input_series = "abycz", $input_charge = "123", $input_decay = "aw")
    {
        $this->m_int = 17999.8;
        $this->m_slp = 9.995216;

        $this->frag_mass_min = 76;
        $this->frag_mass_max = 1800;
        $this->frag_mass_cut = 240;

        $this->setPrecursorMassBin(0.05);

        parent::__construct($input_series, $input_charge, $input_decay);
    }

    public function setMassRange($min, $max)
    {
        if ($min < $max) {
            $this->frag_mass_min = $min;
            $this->frag_mass_max = $max;
        }
    }

    public function setMassCutoff($value = NULL)
    {
        if (! is_nan($value))
            $this->frag_mass_cut = $value;
    }

    public function getFragmentWordArray($aa)
    {
        $array = $this->getFragmentMassArray($aa);

        /*
         * remove values with highly variable mass defect differences
         */
        $array = array_filter($array, function ($x) {
            return ($x > $this->frag_mass_min & $x < $this->frag_mass_max);
        });

        return array_map(array(
            $this,
            'massToHash'
        ), $array);
    }

    public function massToIndex(float $float, $z_detect = false)
    {
        if ($this->frag_mass_cut >= $float)
            return round(($float * $this->m_slp + $this->m_int) / 5) * 5;

        $int = round(($float * $this->m_slp + $this->m_int) / 5) * 5 / 10;

        if (is_true($z_detect) and ! is_int($int))
            $float = chargeMass(neutralMass($float, 2), 1);

        return round(($float * $this->m_slp + $this->m_int) / 5) * 5;
    }

    public function massToHash(float $float, $z_detect = false)
    {
        $int = $this->massToIndex($float, $z_detect = false);
        return base26_encode($int);
    }

    public function precursorMassToHash(float $float)
    {
        $int = round($float * (1 / $this->precursor_mass_bin)) + 18000;
        return 'x' . base26_encode($int);
    }

    public function precursorMassSpread(float $mass, float $span = 0.1)
    {
        $int = max(3, ceil($span / $this->precursor_mass_bin));
        if ($int % 2 == 0)
            $int ++;

        $range = range(0, $int - 1);
        $median = array_median($range);

        foreach ($range as $n => $value) {
            $range[$n] = ($value - $median) * $this->precursor_mass_bin + $mass;
        }
        return $range;
    }

    public function setPrecursorMassBin(float $float)
    {
        $this->precursor_mass_bin = $float;
    }

    public function hashToMass($hash)
    {
        $int = base26_decode($hash);
        return ($int - $this->m_int) / $this->m_slp;
    }

    public function hashOverlap(string $peptide, array $masses)
    {

        /*
         * create the mass hash array
         */
        $mass_hash = array_values(array_map(array(
            $this,
            'massToHash'
        ), ($masses)));

        $mass_array = array_combine($mass_hash, $masses);

        /*
         * create the peptide hash array
         */
        $pept_masses = $this->getFragmentMassArray($peptide);

        $pept_hash = array_values(array_map(array(
            $this,
            'massToHash'
        ), ($pept_masses)));

        $pept_hash_l = sizeof($pept_hash);
        foreach ($pept_hash as $n => $val) {
            $pept_wght[$val] = $n / $pept_hash_l;
        }
        $pept_array = array_combine($pept_hash, $pept_masses);

        /*
         * compute the overlap
         */
        $overlap = array_intersect_key($mass_array, $pept_array);
        foreach ($overlap as $key => $value) {
            $overlap[$key] = abs($value - $pept_array[$key]);
        }

        $overlaps = array_intersect_key($pept_wght, $mass_array);

        $out = [
            'overlap' => array_values($overlap),
            'overlap_score' => array_sum($overlaps)
        ];
        return $out;
    }
}

?>