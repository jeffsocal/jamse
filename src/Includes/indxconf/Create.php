<?php
$exe->cli->header("CREATE WORKING DIRS");

$contents = file_get_contents(get_include_path() . 'dat/jvlnse/peptides.conf');
$contents = preg_replace("/\#INDEX_NAME\#/", $index_name, $contents);

$wrt->writeFile('./conf/i_peptides.conf', $contents);
$exe->cli->message("write", "peptide indexer config");

// /*
// * write out the shell start/stop
// */
// $shell_start = "#!/bin/bash" . PHP_EOL . PHP_EOL;
// $shell_start .= "eval searchd --config peptides.conf" . PHP_EOL;
// $wrt->writeFile('./conf/start_peptides.sh', $shell_start);
// $exe->cli->message("write", "start searchd");
// $wrt->writeFile('./conf/stop_peptides.sh', preg_replace("/searchd/", "searchd --stop", $shell_start));
// $exe->cli->message("write", "stop searchd");

/*
 * CREATE INDEXER BASH
 */
$dirs = [
    './log',
    './data',
    './indx'
];
foreach ($dirs as $dir) {
    if (file_exists($dir))
        continue;
    $exe->cli->message("mkdir", $dir);
    mkdir($dir);
}

?>