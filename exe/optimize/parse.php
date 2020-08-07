#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
searches the jaevalen engine for matching peptides based on the
precursor and fragments masses";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyTime\ProgressTimer;
use ProxyTime\Timer;
use Jvln\Engines\JvlnQuery;
use Jvln\System\JsonConfig;
use ProxyIO\Cli;
use Jvln\Database\Spectra;
use Jvln\Database\PeptideMatches;
use ProxyIO\File\Delim\Delim;
use ProxyIO\File\Delim\WriteDelim;

$cli = new Cli();
$pgt = new ProgressTimer();
$csv = new Delim();
/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v0.9.190917');
$cli->Help->setUsage('-f <file_id>');
$cli->Help->setExample('-db yeast -nph 5');

$ui_path = $cli->getVar('path', "/var/jvlnse/local/optimize/", 'dir path');
/*
 * display the help if needed
 */
$cli->Help->displayHelp();

$files = list_files($ui_path);

foreach ($files as $file) {

    $pgt->cli->message("read file", $file);
    $dat = $csv->read($ui_path . $file);
    $dat['bm15'] = [];
    $dat['bm25a'] = [];
    $dat['field_mask'] = [];
    $dat['doc_word_count'] = [];

    $pgt->start('decode json', table_length($dat));
    foreach ($dat['factors'] as $i => $v) {
        $pgt->barPercent();
        $data = json_decode(str_replace("'", '"', $v), true);

        $dat['bm15'][$i] = $data['bm15'];
        $dat['bm25a'][$i] = $data['bm25a'];
        $dat['field_mask'][$i] = $data['field_mask'];
        $dat['doc_word_count'][$i] = $data['doc_word_count'];
        foreach ($data['fields'][0] as $name => $val) {

            if (! key_exists($name, $dat))
                $dat[$name] = [];

            $dat[$name][$i] = $val;
        }
        // unset($data['words']);
        // print_r($data);
        // if ($i == 5)
        // exit();
    }

    unset($dat['factors']);
    $dat = table_fillCol($dat, 'file_id', str_replace(".csv", "", $file));

    $csv->writeCsvFile($ui_path . "parsed/" . $file, $dat);
}

?>