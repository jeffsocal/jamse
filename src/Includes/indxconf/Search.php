<?php
$exe->cli->header("BUILD SEARCH CONF");

$search_main = file_get_contents(get_include_path() . 'dat/jvlnse/search_apex.conf');
$search_main = preg_replace("/\#INDEX_NAME\#/", $index_name, $search_main);

/*
 * list sharded files
 */
$index_sets = [
    'peptides',
    'homology',
    'fragments'
];

$header = str_repeat("#", 77) . "\n";

$index_source = '';
foreach ($index_sets as $index_set) {
    $index_source .= $header . "# " . $index_set . "\n" . $header;
    $index_mapped = "index " . $index_set . "\n{\n\ttype\t= distributed\n";

    $files = list_files('./indx/', $index_set . ".*\.spd$", true);

    if ($index_set == 'fragments' & sizeof($svr_ns) > 1) {
        foreach ($svr_ns as $svr_n) {
            $index_mapped .= "\tagent_persistent\t= dist" . $svr_n . ":" . $index_port . ":fragments\n";
        }
    } else {
        foreach ($files as $file) {
            $file = str_replace('.spd', '', $file);
            $file = preg_replace('/^\/.*\/jvlnse/', '/var/jvlnse', $file);
            $index_source .= "index " . basename($file) . "\n{\n\tpath\t= " . $file ;
            $index_source .= "\n\tmlock\t= 1\n}\n\n";

            // $index_mapped .= "\tagent\t= localhost:" . $index_port . ":" . basename($file) . "\n";
            $index_mapped .= "\tlocal\t= " . basename($file) . "\n";
        }
    }
    if (count($files) > 1)
        $index_source .= $index_mapped . "}\n";
}
$search_main = preg_replace("/\#NTHREADS\#/", $ui_n_threads, $search_main);
$search_main = preg_replace("/\#INDEX_SOURCE\#/", $index_source, $search_main);
$search_main = preg_replace("/\#LISTEN_PORT_API\#/", $index_port, $search_main);
$search_main = preg_replace("/\#LISTEN_PORT_CLI\#/", $index_port + 1, $search_main);

$wrt->writeFile('./conf/searchd_apex_' . $index_name . '_' . $index_port . '.conf', $search_main);

$shell_run = "#!/bin/bash\n#\n";
$shell_run .= "searchd --config " . './conf/searchd_apex_' . $index_name . '_' . $index_port . '.conf';
$wrt->writeFile("./start", $shell_run);
$wrt->writeFile("./stop", $shell_run . " --stop");

// index i_fpep_hsapiens
// {
// type = distributed
// agent = localhost:9309:i_xaa_fpep_ms
// agent = localhost:9311:i_xab_fpep_ms

if (sizeof($svr_ns) > 1) {

    foreach ($svr_ns as $svr_n) {
        $search_serv = file_get_contents(get_include_path() . 'dat/jvlnse/search_dist.conf');
        $search_serv = preg_replace("/\#INDEX_NAME\#/", $index_name, $search_serv);

        $header = str_repeat("#", 77) . "\n";
        $index_set = 'fragments';

        $index_source = '';
        $files = list_files('./indx/', 'svr_' . $svr_n . '.*' . $index_set . ".*\.spd$", true);
        $index_source .= $header . "# " . $index_set . "\n" . $header;
        $index_mapped = "index " . $index_set . "\n{\n\ttype\t= distributed\n";
        foreach ($files as $file) {
            $file = str_replace('.spd', '', $file);
            $index_source .= "index " . basename($file) . "\n{\n\tpath\t= " . $file . "\n}\n\n";

            // $index_mapped .= "\tagent\t= localhost:" . $index_port . ":" . basename($file) . "\n";
            $index_mapped .= "\tlocal\t= " . basename($file) . "\n";
        }
        if (count($files) > 1)
            $index_source .= $index_mapped . "}\n";
        $search_serv = preg_replace("/\#NTHREADS\#/", sizeof($files), $search_serv);
        $search_serv = preg_replace("/\#INDEX_SOURCE\#/", $index_source, $search_serv);
        $search_serv = preg_replace("/\#LISTEN_PORT_API\#/", $index_port, $search_serv);
        $search_serv = preg_replace("/\#LISTEN_PORT_CLI\#/", $index_port + 1, $search_serv);

        $wrt->writeFile('./conf/searchd_svr' . $svr_n . '_' . $index_name . '_' . $index_port . '.conf', $search_serv);
    }

    $clis = [];

    $clis = [
        'mkdir -p dist/apex/indx/',
        'mkdir -p dist/apex/conf/',
        'mv indx/peptides* dist/apex/indx/',
        'mv conf/searchd_apex* dist/apex/conf/'
    ];
    foreach ($svr_ns as $svr_n) {
        $clis = array_merge($clis, [
            'mkdir -p dist/svr' . $svr_n . '/indx/',
            'mkdir -p dist/svr' . $svr_n . '/conf/',
            'mv indx/svr_' . $svr_n . '* dist/svr' . $svr_n . '/indx/',
            'mv conf/searchd_svr' . $svr_n . '* dist/svr' . $svr_n . '/conf/'
        ]);
    }

    $clis = array_merge($clis, [
        'tar -czvf dist_' . $index_name . '.tar.gz ./dist/'
    ]);

    foreach ($clis as $cli) {
        if (is_false($exe->dohup($cli, TRUE)))
            $exe->cli->header("ERROR - EARLY TERMINATION", TRUE);
    }
}

?>