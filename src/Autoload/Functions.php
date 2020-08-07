<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
function include_functions($config = "ini/config.ini")
{
    
    $include_path = get_include_path();
    $ini = parse_ini_file($include_path . $config, TRUE);
    
    $set = array_values($ini['includes']);
    
    if (is_null($set))
        return false;
    
    if (! is_array($set)) {
        $set = array(
            $set
        );
    }
    foreach ($set as $s) {
        $path = $s;
        $ri_dir = new RecursiveDirectoryIterator($include_path . $path);
        $ri_itr = new RecursiveIteratorIterator($ri_dir, RecursiveIteratorIterator::SELF_FIRST);
        $ri_rgx = new RegexIterator($ri_itr, '/\.php$/i', RecursiveRegexIterator::GET_MATCH);
        
        foreach ($ri_rgx as $name => $object) {
            include $name;
        }
    }
}

/*
 * envoke the includes function
 */

include_functions();


