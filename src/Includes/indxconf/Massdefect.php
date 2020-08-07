<?php
$exe->cli->header("MASS DEFECT MODELING");

$files = scandir('./data/');
$files = array_filter($files, function ($x) {
    return preg_match("/^x.+_peptides$/", $x);
});

if (count($files) > 0) {
    $pids = [];
    foreach ($files as $file) {
        $cli = 'mdefect -f data/' . $file . ' > data/' . $file . '_mdfct';
        if (is_false($exe->nohup($cli, TRUE)))
            $exe->cli->header("ERROR - EARLY TERMINATION", TRUE);
        
        $pids[] = $exe->getPID();
    }
    
    $exe->monitorPIDs($pids, 'parallel mdefect');
} else {
    $exe->cli->header("ERROR - NO PEPTIDE FILES", TRUE);
}

/*
 * dump dedupe peptides
 */
$clis = [
    'cat data/x*peptides*mdfct > data/peptides-mdfct',
    'rm data/x*peptides*mdfct'
];

foreach ($clis as $cli) {
    if (is_false($exe->dohup($cli, TRUE)))
        $exe->cli->header("ERROR - EARLY TERMINATION", TRUE);
}

?>