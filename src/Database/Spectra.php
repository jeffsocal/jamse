<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace Jvln\Database;

class Spectra extends Methods
{

    public function __construct(string $server, int $port = 3306)
    {
        parent::__construct($server, $port);
    }

    public function putFileData(array $json)
    {
        if (strlen($json['file']) > 255)
            return FALSE;

        if (strlen($json['file_id']) > 45)
            return FALSE;

        if (strlen($json['path']) > 255)
            return FALSE;

        if (! is_int($json['scans_n']))
            return FALSE;

            
        /*
         * put the file data into the database
         */
        $this->Prepared->setStatement("INSERT INTO jvln.files 
                             (file_name, file_id, path, scans_n) 
                             VALUES (?,?,?,?)");

        $this->Prepared->addVarible('file_name', $json['file'], 's');
        $this->Prepared->addVarible('file_id', $json['file_id'], 's');
        $this->Prepared->addVarible('path', $json['path'], 's');
        $this->Prepared->addVarible('scans_n', $json['scans_n'], 'i');

        $this->Prepared->paramPut();

        $file_pk = $this->Prepared->getLastInsertID();
        
        /*
         * put the file data into the database
         */
        $this->Prepared->setStatement("INSERT INTO jvln.experiments
                             (file_pk, experiment_name)
                             VALUES (?,?)");
        
        $this->Prepared->addVarible('file_pk', $file_pk, 's');
        $this->Prepared->addVarible('experiment_name', $json['experiment_name'], 's');
        
        $this->Prepared->paramPut();
        
        return $file_pk;
    }

    function putScanData(array $table, array $method, string $file_pk)
    {
        $split = 500;
        $method_pk = $this->putMethods($method);

        if (is_false($method_pk))
            return FALSE;

        $tl = table_length($table);
        $table['file_pk'] = array_fill(0, $tl, $file_pk);
        $table['method_pk'] = array_fill(0, $tl, $method_pk);

        $table['peaks'] = array_map('json_encode', $table['peaks']);
        $table['spectrum'] = array_map('json_encode', $table['spectrum']);

        $this->Simple->setTable('scans');

        $head = table_header($table);
        $table_index_chunks = array_chunk(array_keys($table[$head[0]]), $split);

        foreach ($table_index_chunks as $table_indexs) {

            $sql_string = "INSERT INTO jvln.scans\n";
            $sql_string .= "(" . preg_replace("/\"{2}/", "NULL", array_toString($head)) . ")\n";
            $sql_string = str_replace('"', '', $sql_string);
            $sql_string .= " VALUES\n";

            foreach ($table_indexs as $i) {
                $sql_string_line = "(";
                foreach ($head as $name) {
                    if (preg_match("/^(peaks|spectrum)$/", $name)) {
                        $val = 'COMPRESS(\'' . $table[$name][$i] . '\')';
                    } else {
                        $val = "NULL";
                        if ($table[$name][$i] != '')
                            $val = '"' . $table[$name][$i] . '"';
                    }

                    $sql_string_line .= $val . ",";
                }
                $sql_string_line = trim($sql_string_line, ",");
                $sql_string .= $sql_string_line . "),\n";
            }
            $sql_string = trim($sql_string, ",\n") . ";";

            $return = $this->Simple->sqlPut($sql_string);

            if (is_false($return)) {
                return false;
            }

            // sleep 2 sec to rest the mysql processes
            // sleep(1);
        }

        return true;
    }

    function getScans(string $file_id, int $limit = 0)
    {
        $query = "SELECT scans.pk scan_pk, scan_id, 
                  rt_sec, pre_mz, pre_z, 
                  UNCOMPRESS(peaks) peaks
                  FROM jvln.files, jvln.scans
                  WHERE jvln.files.pk = jvln.scans.file_pk
                  AND n_peaks >= 4
                  AND files.file_id = ?";

        if ($limit > 0)
            $query .= " LIMIT ?";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('file_id', $file_id);

        if ($limit > 0)
            $this->Prepared->addVarible('limit', $limit);

        return $this->Prepared->paramGet();
    }

    function getScan(string $file_id, int $scan_id = 1)
    {
        $query = "SELECT scans.pk scan_pk, scan_id, 
                  rt_sec, pre_mz, pre_z, 
                  UNCOMPRESS(peaks) peaks
                  FROM jvln.files, jvln.scans
                  WHERE jvln.files.pk = jvln.scans.file_pk
                  AND files.file_id = ?
                  AND scans.scan_id = ?";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('file_id', $file_id);
        $this->Prepared->addVarible('scan_id', $scan_id);

        return $this->Prepared->paramGet();
    }
}

?>