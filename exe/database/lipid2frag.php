#!/usr/bin/env php
<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * pep2ptm - apply ptms to digested peptides
 *
 * SYNOPSIS
 * pep2ptm -i <index_name> -f <ptm_file>
 *
 * DESCRIPTION
 *
 * COMMAND LINE
 *
 * -f : <full_path_to_file>
 *
 * -i : <path_to_ini>
 *
 * EXAMPLE
 *
 * ptm2frag -f /var/sphinx/yeast/data/xaa_peptides_ptm_acccd
 *
 */
use ProxyIO\Cli;
use ProxySci\Omics\Proteomics\Peptide;
use Jvln\Engines\TandemHash;

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

$cli = new Cli();
$ui_file = $cli->getVar('f');

// $ui_path_json = list_files('./', ".json", TRUE);
// if (count($ui_path_json) == 0)
// $cli->message('no config file found', "", TRUE);
// $ui_path_json = $ui_path_json[0];
$json_out = [];
$mz_frag = [];
$handle = fopen($ui_file, "r");
if ($handle) {
    while (($line = fgets($handle)) !== false) {
        if (preg_match("/^\d+\.\d+/", $line)) {
            $mz_frag[] = trim(preg_replace("/\s+.+/", '', $line));
        }

        if (preg_match("/^MW/", $line)) {
            $mz_prec = trim(preg_replace("/^[A-Z\:]+\s+/", '', $line));
        }
        if (preg_match("/^Name/", $line)) {
            $name = trim(preg_replace("/^[a-zA-Z\:]+\s+/", '', $line));
            $names = explode(";", $name);
            unset($names[2]);

            $name = array_tostring($names, '', '');
        }
        if (preg_match("/^Comment/", $line)) {
            $comment = trim(preg_replace("/^[a-zA-Z\:]+\s+/", '', $line));
        }
        if (preg_match("/^\s*$/", $line)) {
            $frag_hash = md5(array_tostring($mz_frag, ' ', ''));

            $json_out[$frag_hash] = [
                'fragments' => $mz_frag,
                'mass' => $mz_prec
            ];

            $json_out[$frag_hash]['molecule'][] = [
                'name' => $name,
                'comment' => $comment
            ];

            $mz_frag = [];
        }
    }

    fclose($handle);
} else {
    echo "FAILED TO READ LINE";
}
echo count($json_out) . PHP_EOL;

file_put_contents($ui_file . ".json", json_encode(array_values($json_out), JSON_NUMERIC_CHECK));

?>
