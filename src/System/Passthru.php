<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace Jvln\System;

use ProxyIO\File\Log;
use ProxyTime\Timer;
use ProxyIO\cli;

class Passthru extends Log
{

    private $stdout;

    private $stderr;

    private $pid;

    private $sid;

    public $cli;

    function __construct()
    {
        parent::__construct('cmd');
        $this->cli = new Cli();
    }

    function nohup($command, $print = TRUE)
    {
        /*
         * send execution to backgroud
         */
        return $this->run($command, $print);
    }

    function dohup($command, $print = TRUE)
    {
        /*
         * monitor execution
         */
        if (is_false($this->run($command, $print)))
            return FALSE;

        $this->monitorPIDs($this->getPID(), $command, $print);

        return TRUE;
    }

    function fedup($command, $print = FALSE)
    {
        $tmr = new Timer();
        $this->stdout = [];
        $this->stderr = 0;
        $this->sid = uniqid();
        /*
         * monitor execution
         */
        $this->addToLog(__METHOD__, $command);
        $this->sid = uniqid();

        if (is_true($print))
            $this->cli->flush('eid: ' . $this->sid, substr(date("H:s:m") . " " . $command, 0, $this->cli->getScopeLength() - $this->cli->getMessageLength() - 1));

        exec($command, $this->stdout, $this->stderr);

        $this->stdout[] = $tmr->timeinstr(TRUE);

        if ($this->stderr != 0 || preg_grep("/FAIL/i", $this->stdout)) {
            $this->cli->flusheol("FAILED: " . $tmr->timeinstr(), substr("/tmp/" . $this->sid . " < " . $command, 0, $this->cli->getScopeLength() - $this->cli->getMessageLength() - 1));
            file_put_contents("/tmp/" . $this->sid, array_tostring($this->stdout, "\n", ""));
            return FALSE;
        } else {
            $this->cli->flusheol("SUCCESS: " . $tmr->timeinstr(), substr($command, 0, $this->cli->getScopeLength() - $this->cli->getMessageLength() - 1));
            return TRUE;
        }
    }

    private function run($command, $print = TRUE)
    {
        $this->stdout = [];
        $this->stderr = 0;
        $this->sid = uniqid();

        /*
         * all cli are nohup'd
         */
        $command = 'nohup ' . $command;
        /*
         * if the stdout is not destine for file the capture in tmp
         */
        if (! preg_match("/\s\>\s/", $command))
            $command .= ' > /dev/null';
            
            $command .= ' 2>&1 & echo $!;';

        $this->addToLog(__METHOD__, $command);

        /*
         * run command and capture standard out
         */
        exec($command, $this->stdout, $this->stderr);

        $this->pid = $this->stdout[0];

        // if(!is_false($nohup) & $this->stderr != 0)

        if ($this->stderr != 0)
            return FALSE;
        else
            return TRUE;
    }

    function getStdout()
    {
        return $this->stdout;
    }

    function getPID()
    {
        return $this->pid;
    }

    function getSID()
    {
        return $this->sid;
    }

    function monitorPIDs($pids, string $string = 'exe', bool $print = TRUE)
    {
        if (! is_array($pids))
            $pids = array(
                $pids
            );

        $time_now = $time_start = date('U');
        $n_pids = count($pids);

        $this->cli->flush('exe: starting ' . $n_pids, $string);

        $stdout = [];
        $pids = array_tostring($pids, ' ', '');
        while (2 > 1) {
            $time_now = date('U');
            $elapsed = $time_now - $time_start;

            $stdout = [];
            $sterr = [];
            exec('ps -ylp ' . $pids, $stdout, $sterr);

            if (is_true($print)) {
                $str_left = 'exe: ' . max(0, (count($stdout) - 1)) . '/' . $n_pids;
                $str_rght = '  ' . time_toString($elapsed);
                $this->cli->flush($str_left . $str_rght, $string);
            }

            if (count($stdout) < 2)
                break;

            // $header = preg_split("/\s+/", $stdout[0]);
            // $table = [];
            // unset($stdout[0]);
            // foreach ($stdout as $n => $row) {
            // $row = preg_split("/\s+/", $row);
            // $rowl = array_slice($row, 0, count($header) - 1);
            // $rowl[] = array_tostring(array_slice($row, count($header) - 1), ' ', '');
            // $table[$n] = array_combine($header, $rowl);
            // }

            sleep(1);
        }

        if (is_true($print))
            echo PHP_EOL;

        $this->cli->forceLineCount();

        return TRUE;
    }
}
?>