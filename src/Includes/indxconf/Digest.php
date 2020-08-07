<?php
/*
 * digest protein, index peptides, dump dedupe peptides
 */
$exe->cli->header("DIGEST PROTEINS");
$cmds = [
    'indexer --config conf/i_peptides.conf i_peptides_' . $index_name,
    'indextool --config conf/i_peptides.conf --dumpdict i_peptides_' . $index_name . ' > ./data/peptides',
    'tail -n +8 ./data/peptides > ./data/peptides_tmp',
    'rm ./data/peptides',
    'mv ./data/peptides_tmp  ./data/peptides'
];

foreach ($cmds as $cmd) {
    if (is_false($exe->dohup($cmd, TRUE)))
        $exe->cli->header("ERROR - EARLY TERMINATION", TRUE);
}

?>