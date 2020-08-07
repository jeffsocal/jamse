#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
generates the sphinx database given config json file";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyIO\Cli;

$cli = new Cli();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('RELEASE v1.2.19');
$cli->Help->setUsage('');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

function human_filesize($bytes, $decimals = 2)
{
    $size = array(
        'B',
        'kB',
        'MB',
        'GB',
        'TB',
        'PB',
        'EB',
        'ZB',
        'YB'
    );
    $factor = floor((strlen($bytes) - 1) / 3);
    return sprintf("%.{$decimals}f", $bytes / pow(1024, $factor)) . @$size[$factor];
}

function get_filelength($file)
{
    if (! file_exists($file))
        return false;

    exec("wc -l " . $file, $stdout, $stderr);
    return preg_replace("/\s.+/", '', $stdout[0]);
}

function get_fastacount($file)
{
    if (! file_exists($file))
        return false;
    exec('eval cat ' . $file . ' | grep ">" | wc -l', $stdout, $stderr);
    return preg_replace("/\s.+/", '', $stdout[0]);
}

echo PHP_EOL;
$cli->header("JAEVALEN DATABASE STATS");
$cli->header("CONFIG VALUES");

$json_file_path = list_files('./', ".json", TRUE);
if (count($json_file_path) == 0)
    $cli->message('no config file found', "", TRUE);
$json_file_path = $json_file_path[0];

$ini = parse_json_file($json_file_path);
$cli->message('config path', $json_file_path);

/*
 * SET THE PARAMS
 */
$index_name = $ini['index']['name'];
$unimod_set_name = $ini['unimod']['set_name'];
$cli->message('index name', $index_name);

$prot_size = 0;
foreach ($ini['fasta']['files'] as $f) {

    $file_path = $ini['fasta']['path'] . $f;
    $prot_size_this = 0;

    if (file_exists($file_path))
        $prot_size_this = get_fastacount($file_path);

    $cli->message('uniprot | ' . number_format($prot_size_this), $f);
    $prot_size += $prot_size_this;
}

$mods_size = 0;
foreach ($ini['unimod']['modification_set'] as $name => $vals) {
    $mods_size += count($vals['applied_sites']);
    $cli->message('unimod', $name . ' | ' . array_tostring($vals['applied_sites'], ' ', ''));
}

$cli->message('unimod applied', $mods_size);
$cli->message('unimod concurrent', $ini['unimod']['max_concurrent']);

$cli->header("INDEX VALUES");
$cli->message('data/proteins', number_format($prot_size));
$cli->message('data/peptides', number_format(get_filelength('data/peptides')));

/*
 * report the number of peptides w/ ptms
 */
$ptms_size = 1;
foreach (list_files('data/', 'fragments', TRUE) as $f) {
    if (file_exists($f))
        $ptms_size += get_filelength($f);
}
$cli->message('data/peptide-ptms', number_format($ptms_size));

$data_size = 0;
foreach (scandir('./data/') as $f)
    $data_size += filesize('./data/' . $f);

$indx_size = 0;
foreach (scandir('./indx/') as $f)
    $indx_size += filesize('./indx/' . $f);

$cli->message("data/", human_filesize($data_size, 1));
$cli->message("indx/", human_filesize($indx_size, 1));

$cli->header("COMPLETED");
echo PHP_EOL;
?>