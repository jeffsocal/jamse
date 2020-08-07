<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace Jvln\System;

class JsonConfig
{

    protected $file_contents;

    protected $file_path;

    public function __construct($file_path = false)
    {
        $this->file_path = false;
        $this->file_contents = false;
        if (is_false($file_path))
            $file_path = list_files('./', ".json", TRUE);
        else
            $file_path = list_files($file_path, ".json", TRUE);

        if (is_null($file_path))
            return false;

        $this->file_path = $file_path[0];

        if (is_false(file_exists($this->file_path)))
            return false;

        $this->file_contents = json_decode(file_get_contents($this->file_path), true, 10);
    }

    public function getConfigFile()
    {
        return $this->file_path;
    }

    public function getVariable(string $group = null, string $variable = null, $value = null)
    {
        $l1 = $group;
        $l2 = $variable;

        if (is_null($this->file_contents))
            return $value;

        if (is_null($l1))
            return $value;

        if (! key_exists($l1, $this->file_contents))
            return $value;

        if (is_null($l2))
            return $this->file_contents[$l1];

        if (! key_exists($l2, $this->file_contents[$l1]))
            return $value;

        return $this->file_contents[$l1][$l2];
    }
}
?>
