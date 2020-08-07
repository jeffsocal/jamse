<?php
$exe->cli->message("index fragmentation");

$shard_main = file_get_contents(get_include_path() . 'dat/jvlnse/shard_main.conf');

$shard_index = file_get_contents(get_include_path() . 'dat/jvlnse/shard_index.conf');
$shard_index = preg_replace("/\#INDEX_NAME\#/", $index_name, $shard_index);
$shard_index = preg_replace("/\#JVLN_EXE\#/", 'ptm2frag', $shard_index);
$shard_index = preg_replace("/\#MIN_WORD_LENGTH\#/", '3', $shard_index);

$sever_prefix = 'svr_' . $svr_n;
$shard_prefix = $sever_prefix . '_' . $unimod_set_name;

/*
 * shard peptides-ptms
 */
$clis = [
    'split -n r/' . $ui_n_threads . ' data/' . $shard_prefix . ' ' . $shard_prefix . '_ --additional-suffix=_fragments',
    'mv ' . $shard_prefix . '*fragments data/',
    'rm indx/' . $shard_prefix . '*fragments'
];

foreach ($clis as $cli) {
    if (is_false($exe->dohup($cli, TRUE)))
        $exe->cli->header("ERROR - EARLY TERMINATION", TRUE);
}

/*
 * list sharded files
 */
$files = list_files('./data/', $shard_prefix . ".*fragments$");

$clis = [];
if (count($files) > 0) {
    foreach ($files as $file) {

        $this_shard_index = $shard_index;

        $this_shard_index = preg_replace("/\#SHARD_NAME\#/", $file, $this_shard_index);
        $this_shard_index = preg_replace("/\#FILE_PATH\#/", "data/" . $file, $this_shard_index);

        $shard_main = preg_replace("/\#INDEX_SOURCE#/", $this_shard_index, $shard_main);

        $clis[] = 'indexer --config conf/' . $sever_prefix . '_fragments.conf i_' . $file;
    }

    /*
     * write out the master config file
     */
    $wrt->writeFile('./conf/' . $sever_prefix . '_fragments.conf', $shard_main);

    /*
     * run the indexer on the peptide-ptm file
     */
    $pids = [];
    foreach ($clis as $cli) {
        if (is_false($exe->nohup($cli, TRUE)))
            $exe->cli->header("ERROR - EARLY TERMINATION", TRUE);
        $pids[] = $exe->getPID();
    }
    $exe->monitorPIDs($pids, 'exe || fragments');
} else {
    $exe->cli->header("ERROR - NO PEPTIDE-PTM FILES", TRUE);
}

?>
