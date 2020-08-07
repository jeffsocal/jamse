<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace Jvln\Engines;

use ProxyIO\File\Log;
use mysqli;

class Spectrum
{

    private $query_peaks;

    private $query_premz;

    private $query_prez;

    private $query_preia;

    private $query_limit;

    private $jvln_s;

    private $jvln_p;

    private $jvln_db;

    private $jvln_un;

    private $jvln_pw;

    private $jvln_ranker;

    protected $index;

    public function __construct($server = '127.0.0.1', $port = 9307)
    {
        $this->jvln_s = $server;
        $this->jvln_p = $port;
        $this->jvln_db = NULL;
        $this->jvln_un = 'root';
        $this->jvln_pw = '';
        $this->jvln_ranker = "bm25*sum(word_count)";
    }

    function setIndex(string $index)
    {
        $this->index = $index;
    }

    function setRankerFunction(string $function)
    {
        $this->jvln_ranker = $function;
    }

    function search(array $peaks, float $pre_mz, $pre_zs = NULL, $pre_ia = NULL, int $limit = 50)
    {
        $this->query_peaks = $peaks;
        $this->query_premz = $pre_mz;
        $this->query_prez = $pre_zs;
        $this->query_preia = $pre_ia;
        $this->query_limit = $limit;

        $hits = $this->queryDB();

        if (is_false($hits))
            return FALSE;

        return $this->queryMunge($hits);
    }

    function queryDB()
    {
        $log = new Log("sphinx");

        $jvln_lh = $this->query_limit;

        if (count($this->query_peaks) < 2)
            return false;

        if (is_null($this->query_prez) | $this->query_prez == 0)
            $query_pzs = range(1, 3);

        if (! is_array($this->query_prez))
            $query_pzs = [
                $this->query_prez
            ];

        /*
         *
         * SQL SERVER CONN
         *
         */

        // $pgt->message("jvln", "connect");
        $jvln_i = @new mysqli($this->jvln_s, $this->jvln_un, $this->jvln_pw, $this->jvln_db, $this->jvln_p);

        $jvln_i_id = $jvln_i->thread_id;

        if ($jvln_i->connect_errno) {
            $log->addToLog("connect error: " . $jvln_i->connect_errno, $jvln_i->connect_error);
            $log->addToLog("thread id", $jvln_i_id);
        }

        $jvln_data = [];
        /*
         * FOREACH SPECTRA -- START
         */

        // POW(10, (-1 * WEIGHT() / 1000)) score
        $jvln_sq = getSphinxQuery($this->query_peaks, $this->query_premz, $query_pzs, NULL, 0, 0.5);

        $jvln_q = 'SELECT peptide, 
                    WEIGHT() / 10000 as score ';

        // if ($this->jvln_ranker == '1')
        $jvln_q .= ', PACKEDFACTORS({json=1, no_atc=1}) as factors  ';

        $jvln_q .= 'FROM ' . $this->index . "
                    WHERE MATCH('" . $jvln_sq . "')
                    LIMIT " . $jvln_lh . "
                    OPTION ranker=expr('" . $this->jvln_ranker . "')";

        if ($jvln_r = $jvln_i->query($jvln_q)) {
            $jvln_data = $jvln_r->fetch_all(MYSQLI_ASSOC);
            $jvln_r->free();
        }

        /*
         * FOREACH SPECTRA -- END
         */

        $jvln_err = $jvln_i->error;
        $jvln_i->close();
        return $jvln_data;
    }
}

?>