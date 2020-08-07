#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
generates the sphinx database given config json file
databases need to be indexed with a unique 3 digit value.
This will act as the DB port designation in the 55000 range
and allow for multiple DB's to be active on a single machine
since most require minimal RAM.";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyIO\Cli;
use Jvln\System\Passthru;
use Jvln\System\JsonConfig;
use ProxyIO\File\Write;

$cli = new Cli();
$exe = new Passthru();
$jsc = new JsonConfig();
$wrt = new Write();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('RELEASE v1.1.190823');
$cli->Help->setUsage('-m <methods>');
$cli->Help->setExample('-m cre+dig');
$cli->Help->setExample('-m mod+frag+hom');

/*
 * load the variables
 * flag : default-value : description
 */
$ui_n_threads = $cli->getVar('t', 2, 'n threads');
$ui_n_servers = $cli->getVar('s', 1, 'n servers');
$ui_methods = $cli->getVar('m', 'create+digest+modification+fragment+homology', '
<methods> default = create+digest+modification+fragment
-- this notation strings together computational modes
-- with partial spelling accepted (i.e. cre+dig+mod+frag)
-- 
-- create ......... establish the working directories
-- digest ......... protein digestion to pep (w/ indexing)
-- homology ....... create an index of pep for homology searching
-- modification ... peptide expandion by modification
-- fragment ....... peptide-ptm fragmentation (w/ indexing)
-- massdefect ..... calculate mass values for modeling
-- searchd ........ configure the searchd parameters
-- 
--');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

echo PHP_EOL;
$exe->cli->header("CREATE JAEVALEN DATABASE");
$exe->cli->message('ui methods', $ui_methods);
$exe->cli->header("CONFIG VALUES");

/*
 * get INI definitions
 */
if ($jsc->getConfigFile() == false)
    $cli->message('no config file found', "", TRUE);

/*
 * SET THE PARAMS
 */
$index_name = $jsc->getVariable('index', 'name');
$index_port = $jsc->getVariable('index', 'port');
$unimod_modification_set = $jsc->getVariable('unimod', 'modification_set');
$unimod_set_name = $jsc->getVariable('unimod', 'set_name');
$fasta_files = $jsc->getVariable('fasta', 'files');

/*
 * message the user
 */
$exe->cli->message('index name', $index_name);
$exe->cli->message('index port', $index_port);
$exe->cli->message('max indexing threads', $ui_n_threads);
$exe->cli->message('uniprot', trim(array_tostring($fasta_files, ',', ''), '"'));
$exe->cli->message('unimod', array_tostring(array_keys($unimod_modification_set), ',', ''));

if (preg_match("/cre[ate]*/", $ui_methods)) {
    require_once get_include_path() . 'src/Includes/indxconf/Create.php';
}

if (preg_match("/dig[est]*/", $ui_methods)) {
    require_once get_include_path() . 'src/Includes/indxconf/Digest.php';
}

/*
 * propigate servers
 */
$svr_ns = [];
for ($svr = 0; $svr < $ui_n_servers; $svr ++) {
    $svr_ns[] = str_pad($svr, 2, 0, STR_PAD_LEFT);
}

/*
 * determine the optimum sharding -- required for mods|homology|mass-defect
 */
if (preg_match("/mod[ification]*|hom[ology]*|m*(def|mas)[sdefect]*/", $ui_methods)) {

    $exe->cli->header("SHARD DB");
    /*
     * shard by server / then by cores
     */
    $cmds = [
        'rm data/*_peptides',
        'split -dn r/' . $ui_n_servers . ' ./data/peptides svr_ --additional-suffix=_peptides',
        'mv svr*peptides' . ' data/'
    ];
    /*
     * then shard by cores
     */
    foreach ($svr_ns as $svr_n) {
        $cmds = array_merge($cmds, [
            'split -n r/' . $ui_n_threads . ' ./data/svr_' . $svr_n . '_peptides  svr_' . $svr_n . '_peptides_',
            'mv svr_' . $svr_n . '_peptides*' . ' data/',
            'rm data/svr_' . $svr_n . '_peptides'
        ]);
    }
    foreach ($cmds as $cmd) {
        if (is_false($exe->dohup($cmd, true))) {
            $exe->cli->header("ERROR - EARLY TERMINATION", true);
        }
    }
}


foreach ($svr_ns as $svr_n) {

    $exe->cli->header('CONFIG SERVER '. $svr_n);
    /*
     *
     */
    if (preg_match("/hom[ology]*/", $ui_methods)) {
        require get_include_path() . 'src/Includes/indxconf/Homology.php';
    }

    /*
     * server: processes --> recombines
     */
    if (preg_match("/mod[ification]*/", $ui_methods)) {
        require get_include_path() . 'src/Includes/indxconf/Modification.php';
    }

    if (preg_match("/m*(def|mas)[sdefect]*/", $ui_methods)) {
        require get_include_path() . 'src/Includes/indxconf/Massdefect.php';
    }

    /*
     * server: shards --> processes --> recombines
     */
    if (preg_match("/frag[ment]*/", $ui_methods)) {
        require get_include_path() . 'src/Includes/indxconf/Fragment.php';
    }
}

/*
 * clean up shard
 */
$clis = [
    'rm data/svr*'.$unimod_set_name
];
foreach ($clis as $cli) {
    if (is_false($exe->dohup($cli, TRUE)))
        $exe->cli->header("ERROR - EARLY TERMINATION", TRUE);
}

require_once get_include_path() . 'src/Includes/indxconf/Search.php';

$exe->cli->header("COMPLETED");

echo PHP_EOL;

?>