<?php
/*
 * kick off the propigation of ptms
 */
$exe->cli->message("apply ptms to peptides");

$files = list_files('./data/', "^svr_" . $svr_n . "_peptides_[a-z]{2}$");

/*
 * run files in parallel
 */
if (count($files) > 0) {
    $pids = [];
    foreach ($files as $file) {
        $cli = 'pep2ptm -f data/' . $file . ' > data/' . $file . '_ptm';
        if (is_false($exe->nohup($cli, TRUE)))
            $exe->cli->header("ERROR - EARLY TERMINATION", TRUE);

        $pids[] = $exe->getPID();
    }

    $exe->monitorPIDs($pids, 'exe || pep2ptm');
} else {
    $exe->cli->header("ERROR - NO PEPTIDE FILES", TRUE);
}

/*
 * dump dedupe peptides
 */
$clis = [
    'cat data/svr_' . $svr_n . '_peptides*ptm > data/svr_' . $svr_n . '_' . $unimod_set_name,
    'rm data/svr_' . $svr_n . '_peptides*ptm'
];

foreach ($clis as $cli) {
    if (is_false($exe->dohup($cli, TRUE)))
        $exe->cli->header("ERROR - EARLY TERMINATION", TRUE);
}

?>