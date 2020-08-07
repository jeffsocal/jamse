<?php
$exe->cli->message("index homology");

$shard_main = file_get_contents(get_include_path() . 'dat/jvlnse/shard_main.conf');

$shard_index = file_get_contents(get_include_path() . 'dat/jvlnse/shard_index.conf');
$shard_index = preg_replace("/\#INDEX_NAME\#/", $index_name, $shard_index);
$shard_index = preg_replace("/\#JVLN_EXE\#/", 'pep2sap', $shard_index);
$shard_index = preg_replace("/\#MIN_WORD_LENGTH\#/", '3', $shard_index);


$sever_prefix = 'svr_' . $svr_n;
$shard_prefix = $sever_prefix . '_' . $unimod_set_name;

/*
 * list sharded files
 */
$files = list_files('./data/', "^svr_" . $svr_n . "_peptides_[a-z]{2}$");

if (count($files) > 0) {
    $clis = [];
    foreach ($files as $file) {

        $this_shard_index = $shard_index;

        $shard_name = preg_replace("/peptides/", "homology", $file);
        $this_shard_index = preg_replace("/\#SHARD_NAME\#/", $shard_name, $this_shard_index);
        $this_shard_index = preg_replace("/\#FILE_PATH\#/", "data/" . $file, $this_shard_index);

        $shard_main = preg_replace("/\#INDEX_SOURCE#/", $this_shard_index, $shard_main);

        $clis[] = 'indexer --config conf/' . $sever_prefix . '_homology.conf i_' . $shard_name;
    }

    /*
     * write out the master config file
     */
    $wrt->writeFile('./conf/' . $sever_prefix . '_homology.conf', $shard_main);

    /*
     * run the indexer on the peptide-ptm file
     */
    $pids = [];
    foreach ($clis as $cli) {
        if (is_false($exe->nohup($cli, TRUE)))
            $exe->cli->header("ERROR - EARLY TERMINATION", TRUE);
        $pids[] = $exe->getPID();
    }
    $exe->monitorPIDs($pids, 'exe || homology');
} else {
    $exe->cli->header("ERROR - NO PEPTIDE-PTM FILES", TRUE);
}

?>