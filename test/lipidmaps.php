<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
use ProxySci\Omics\Proteomics\Digest;
use ProxySci\Omics\Proteomics\Fasta;
use ProxySci\Omics\Proteomics\Peptide;
use ProxyTime\Timer;
use ProxyIO\File\Read;
use ProxyTime\ProgressTimer;
use ProxyIO\File\Delim\WriteDelim;

$dat = new Read('/data/OpenAccess/LIPID_MAPS/LMSDFDownload12Dec17/LMSDFDownload12Dec17FinalAll.sdf');
$tmr = new Timer();
$dig = new Digest();
$fst = new Fasta();
$amu = new Peptide();
$ptm = new ProgressTimer();
$csv = new WriteDelim();

// echo PHP_EOL;

$data = $dat->getContents();

$keys = [
    'LM_ID',
    'COMMON_NAME',
    'SYSTEMATIC_NAME',
    'SYNONYMS',
    'CATEGORY',
    'MAIN_CLASS',
    'SUB_CLASS',
    'EXACT_MASS',
    'FORMULA',
    'LIPIDBANK_ID',
    'PUBCHEM_SID',
    'PUBCHEM_CID',
    'KEGG_ID',
    'HMDBID',
    'CHEBI_ID',
    'INCHI_KEY',
    'INCHI'
];

$data_chunks = explode("$$$$", $data);

$ptm->proTimerSize(count($data_chunks));
$ptm->proTimerStart('parse data file');
foreach ($data_chunks as $n => $data_chunk) {
    $ptm->proTimerPrint();
    foreach ($keys as $key) {
        $table_data[strtolower($key)][$n] = '';
        if (preg_match_all("/(?<=\>\s\<" . $key . "\>)\n.+/", $data_chunk, $matches) > 0)
            $table_data[strtolower($key)][$n] = trim(str_replace('"', '', $matches[0][0]));
    }
}

print_table($table_data);

$csv->writeCsvFile("/data/OpenAccess/LIPID_MAPS/LMSDF_20171212.csv", $table_data);


?>