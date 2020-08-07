<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace Jvln\Engines;

use Jvln\System\JsonConfig;

class JvlnQuery
{

    private $pre_charge_default;

    private $query_n_must_peaks;

    private $query_quorum;

    private $query_discriminant_fun;

    private $query_return_limit;

    private $query_table;

    private $query_port;

    public function __construct(string $database, string $port)
    {
        $this->query_table = $database;
        $this->query_port = $port;
        $this->setQuorum();
        $this->setNMustHavePeaks();
        $this->setPrecursorChargeDefault();

        // $this->setDiscriminantFunction('(doc_word_count*0.083910471 + bm15*0.003281405 - 2.029638545) * 1000');
        // $this->setDiscriminantFunction('(doc_word_count*6.047 + (doc_word_count/query_word_count)*108.510 + bm15*4.346 + sum(sum_idf)*-1593.438 + sum(wlccs)*-325.969 - 2285.585) * 10');
        /*
         * estimated AUC 0.84 10x10 FCV
         *
         * (Intercept) 10.01712135
         * bm15 -0.01992901
         * tf_idf 6.62916579
         * doc_word_count 0.06004149
         */
        // $this->setDiscriminantFunction('(10.01712135 + bm15*-0.01992901 + sum(tf_idf)*6.62916579 + doc_word_count*0.06004149) * 1000');

        $this->setDiscriminantFunction('(doc_word_count * sum(tf_idf) * 1000)');

        $this->setReturnLimit();
    }

    function setReturnLimit(int $int = 3)
    {
        $this->query_return_limit = $int;
    }

    function setQuorum(float $float = 4)
    {
        $this->query_quorum = $float;
    }

    function setNMustHavePeaks(int $int = 0)
    {
        $this->query_n_must_peaks = $int;
    }

    function setPrecursorChargeDefault(int $int = 2)
    {
        $this->pre_charge_default = $int;
    }

    function setDiscriminantFunction(string $formula)
    {
        $this->query_discriminant_fun = $formula;
    }

    function get_fastacount($file)
    {
        if (! file_exists($file))
            return false;
        exec('eval cat ' . $file . ' | grep ">" | wc -l', $stdout, $stderr);
        return preg_replace("/\s.+/", '', $stdout[0]);
    }

    function getDBStats($server)
    {
        $jsc = new JsonConfig('/var/jvlnse/' . $this->query_table);
        $stats = [];

        /*
         * SET THE PARAMS
         */
        $index_name = $jsc->getVariable('index', 'name');
        $index_port = $jsc->getVariable('index', 'port');

        $fasta_files = $jsc->getVariable('fasta', 'files');
        $fasta_path = $jsc->getVariable('fasta', 'path');
        $prot_size = 0;
        foreach ($fasta_files as $f) {
            $file_path_this = $fasta_path . $f;
            $prot_size_this = 0;
            if (file_exists($file_path_this))
                $prot_size_this = $this->get_fastacount($file_path_this);
            $stats['proteins'][] = [
                'source' => $f,
                'size' => $prot_size_this
            ];
        }
        $stats['modifications'] = $jsc->getVariable('unimod');

        // $jvn = new Simple($server, $index_port);
        // $tbl = $jvn->sqlGet("SELECT count(*) n FROM fragments;");
        // $stats['peptides'] = $tbl['n'];

        return $stats;
    }

    function getQuery(array $peaks, float $pre_mz = NULL, array $pre_z = NULL, array $pre_iso = NULL)
    {
        $this_function = $this->query_discriminant_fun;

        // POW(10, (-1 * WEIGHT() / 1000)) score
        $this_match = $this->getMatchQuery($peaks, $pre_mz, $pre_z, $pre_iso);
        $this_query = "SELECT peptide,
                   WEIGHT() score
                   FROM fragments
                   WHERE MATCH('" . $this_match . "')
                   AND score > 2
                   LIMIT " . $this->query_return_limit . "
                   OPTION ranker=expr('" . $this_function . "')";

        return preg_replace("/\s+/", ' ', $this_query . ";");
    }

    function getPrecursorCount(float $pre_mz = NULL, array $pre_z = NULL, array $pre_iso = NULL)
    {
        $this_function = $this->query_discriminant_fun;

        // POW(10, (-1 * WEIGHT() / 1000)) score
        $this_match = $this->getPrecursorMatch($pre_mz, $pre_z, $pre_iso);
        $this_query = "SELECT count(*) n
                   FROM fragments
                   WHERE MATCH('" . $this_match . "')";

        return preg_replace("/\s+/", ' ', $this_query . ";");
    }

    function getQueryModel(array $peaks, float $pre_mz = NULL, array $pre_z = NULL, array $pre_iso = NULL)
    {
        $this_function = $this->query_discriminant_fun;

        // POW(10, (-1 * WEIGHT() / 1000)) score
        $this_match = $this->getMatchQuery($peaks, $pre_mz, $pre_z, $pre_iso);
        $this_query = "SELECT peptide, WEIGHT() as score,
                   PACKEDFACTORS({json=1}) as factors
                   FROM fragments
                   WHERE MATCH('" . $this_match . "')
                   LIMIT 50
                   OPTION ranker=expr('" . $this->query_discriminant_fun . "')";

        return preg_replace("/\s+/", ' ', $this_query . ";");
    }

    private function getPrecursorMatch(float $pre_mz = NULL, array $pre_z = NULL, array $pre_iso = NULL)
    {
        $frg = new TandemHash();
        /*
         * fill default precursor charge
         */
        if ($pre_z == NULL | $pre_z == 0)
            $pre_z = $this->pre_charge_default;

        if (! is_array($pre_z))
            $pre_z = [
                $pre_z
            ];

        /*
         * fill default precursor isotope
         */
        if ($pre_iso == NULL)
            $pre_iso = 0;

        if (! is_array($pre_iso))
            $pre_iso = [
                $pre_iso
            ];

        $pre_hashes = [];

        if ($pre_mz != NULL) {
            foreach ($pre_z as $z) {
                foreach ($pre_iso as $a) {

                    $pre_bins = $frg->precursorMassSpread(neutralMass($pre_mz, $z) + $a);
                    foreach ($pre_bins as $pre_bin) {
                        $pre_hashes[] = $frg->precursorMassToHash($pre_bin);
                    }
                }
            }
            /*
             * pre cursor mass spread
             */
            $pre_hashes = '"' . array_tostring($pre_hashes, " ", "") . '"/1';
        }
        return $pre_hashes;
    }

    private function getMatchQuery(array $peaks, float $pre_mz = NULL, array $pre_z = NULL, array $pre_iso = NULL)
    {
        $pre_hashes = $this->getPrecursorMatch($pre_mz, $pre_z, $pre_iso);

        /*
         * hash the spectrum peaks
         */
        $frg = new TandemHash();
        $hash_array = array_unique(array_map(array(
            $frg,
            'massToHash'
        ), ($peaks)));

        /*
         * create the must have query
         */
        // $peaks_must = ' ' . array_tostring(array_slice($hash_array, 0, $this->query_n_must_peaks), ' ', '');
        $peaks_must = '';
        if ($this->query_n_must_peaks > 0)
            $peaks_must = ' "' . array_tostring(array_slice($hash_array, 0, $this->query_n_must_peaks), ' ', '') . '"/1';

        /*
         * create the query for the other by quorum
         */
        $hash_array = array_slice($hash_array, $this->query_n_must_peaks);
        // sort($hash_array);

        $peaks_maybe = ' "' . array_tostring($hash_array, ' ', '') . '"/' . $this->query_quorum;

        return $pre_hashes . $peaks_must . $peaks_maybe;
    }
}

?>