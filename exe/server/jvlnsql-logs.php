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

use ProxyTime\Timer;
use ProxyIO\Cli;
use ProxyIO\File\Delim\WriteDelim;

$cli = new Cli();
$tmr = new Timer();
$csv = new WriteDelim();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v0.9.190917');
$cli->Help->setUsage('-f <file_id>');

/*
 * load the variables
 * flag : default-value : description
 */
$file_reg = $cli->getVar('f', NULL, 'file regex');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

if (is_null($file_reg))
    $cli->message("ERROR", "file regex (-f) not set", TRUE);

$files = list_files("./", $file_reg, TRUE);

foreach ($files as $file) {
$table = []; 

    $cont = file_get_contents($file);
    $cont = explode("\n", $cont);
    foreach ($cont as $n => $line) {

        preg_match("/^\[[A-Z][a-z]*.+20[0-9]{2}\]/", $line, $this_date);

        if (sizeof($this_date) == 0)
            continue;

        $this_date = $this_date[0];

        $line = str_replace($this_date, '', $line);

        $this_time = array_slice(explode(' ', $line), 1, 4);

        preg_match_all("/\[[a-z].+?\]/i", $line, $this_refs);

        $table[] = [
            'date' => trim($this_date, "[]"),
            'time_sys' => $this_time[0],
            'time_rel' => $this_time[2],
            'query' => trim($this_refs[0][0], "[]"),
            'table' => trim($this_refs[0][1], "[]")
        ];
    }

    $table = table_invert($table);
    $csv->writeCsvFile(str_replace(".log", "", $file) . '.csv', $table);

}

?>