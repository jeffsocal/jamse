<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace Jvln\Database;

use ProxyIO\Cli;

class Experiments extends Methods
{

    public function __construct(string $server = '127.0.0.1', int $port = 3306)
    {
        parent::__construct($server, $port);
    }

    public function getTableByMatch($str)
    {
        $query = "SELECT
                        experiment_name,
                        file_pk,
                        file_name,
                        file_id,
                        scans_n,
                        CONCAT_WS(' ', experiment_name,
                                file_name,
                                file_id,
                                scans_n) as text
                    FROM
                        jvln.experiments,
                        jvln.files
                    WHERE
                        jvln.experiments.file_pk = jvln.files.pk
                    HAVING
                        text like ?";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('text', "%" . $str . "%");

        $data = $this->Prepared->paramGet();

        if (is_false($data))
            return FALSE;

        return table_dropcols($data, 'text');
    }

    public function deletePSMSByMatch($str)
    {
        $cli = new Cli();
        $table = table_invert($this->getTableByMatch($str));

        foreach ($table as $row) {
            $file_pk = $row['file_pk'];

            $query = "DELETE
                    FROM psms
                    WHERE scan_pk IN 
                    (SELECT pk FROM scans WHERE file_pk = ?);";

            $this->Prepared->setStatement($query);
            $this->Prepared->addVarible('file_pk', $file_pk);

            $status = 'success';
            if ($this->Prepared->paramDel() == FALSE) {
                $status = 'failed';
            }

            $cli->message($status, " deleted:psms " . $row['experiment_name'] . ' ' . $row['file_id']);
        }
    }
}

?>