<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
use ProxyTime\Timer;
use ProxyIO\File\Delim\WriteDelim;

$tmr = new Timer();
$csv = new WriteDelim();

$seq = opts('s');

echo PHP_EOL;
echo $seq . PHP_EOL;

$ini = parse_ini_file(get_include_path() . 'ini/molecular.ini');
$unixml_path = get_include_path() . $ini['unixml_path'];

$xml = file_get_contents($unixml_path);
$p = xml_parser_create();
xml_parse_into_struct($p, $xml, $vals, $index);
xml_parser_free($p);

$array_mods = [
    'Carbamidomethyl',
    'Carbamyl',
    'Acetyl',
    'Hydroxylation',
    'Deamidation',
    'Carboxymethyl',
    'Pyro_glu',
    'Pyro-glu',
    'Phospho',
    'GlyGly',
    'Methyl',
    'Dimethyl'
];

$unimod_mod = [];
foreach ($index['MODIFICATIONS_ROW'] as $ele1 => $wtf) {
    if (key_exists('attributes', $vals[$wtf]))
        $unimod_mod[] = preg_replace("/^\s*|\s*$/", "", $vals[$wtf]['attributes']);
}
$unimod_spe = [];
foreach ($index['SPECIFICITY_ROW'] as $ele1 => $wtf) {
    if (key_exists('attributes', $vals[$wtf]))
        $unimod_spe[] = preg_replace("/^\s*|\s*$/", "", $vals[$wtf]['attributes']);
}

$table_unimod_mod = array_rowtocol($unimod_mod);
$table_unimod_spe = array_rowtocol($unimod_spe);

$table_mod = $table_unimod_mod;
$table_spe = $table_unimod_spe;

$table_out = [];
$table_out['CODE_NAME'] = array_intersect($table_mod['CODE_NAME'], $array_mods);
$table_out['FULL_NAME'] = array_intersect_key($table_mod['FULL_NAME'], $table_out['CODE_NAME']);
$table_out['MONO_MASS'] = array_intersect_key($table_mod['MONO_MASS'], $table_out['CODE_NAME']);
$table_out['RECORD_ID'] = array_intersect_key($table_mod['RECORD_ID'], $table_out['CODE_NAME']);
$table_out['SPE_AMINO'] = [];

foreach ($table_out['RECORD_ID'] as $n => $rec_id) {
    $arr_id = [
        $rec_id
    ];
    $rec_ids = array_intersect($table_spe['MOD_KEY'], $arr_id);
    $spe_class = array_intersect_key($table_spe['CLASSIFICATIONS_KEY'], $rec_ids);
    
    $ptms = array_intersect($spe_class, [2,5,13,6]);
    $spe_amino = array_intersect_key($table_spe['ONE_LETTER'], $ptms);
    
    $table_out['SPE_AMINO'][] = array_tostring(array_unique($spe_amino), ' ', '');
}

$table_out['CODE_NAME'] = array_values($table_out['CODE_NAME']);
$table_out['FULL_NAME'] = array_values($table_out['FULL_NAME']);
$table_out['MONO_MASS'] = array_values($table_out['MONO_MASS']);
$table_out['RECORD_ID'] = array_values($table_out['RECORD_ID']);
$table_out['SPE_AMINO'] = array_values($table_out['SPE_AMINO']);

print_table($table_out);

$csv->writeCsvFile("../dat/unimod/unimod.csv", $table_out);

echo PHP_EOL;
print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);

?>