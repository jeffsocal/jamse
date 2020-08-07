<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace Jvln\System;

use BlakeGardner\MacAddress;
use ProxyIO\Cli;

class Interact extends Cli
{

    // license mac address
    private $jv7s036n;
    
    // license date
    private $jv8nO2Kz;

    function __construct(int $width = 100, string $sep = " ")
    {
        parent::__construct($width, $sep);
        $this->setLicenseHash('WyJERToxNTo1RjpGNzo4MjpGQiIsMTcwMDU4NTIwMCwzMjFd');
//         $this->setLicenseHash('WyJERToxNTo1RjpGNzo4MjpIwMCwzMjFd');
        $this->killApplicaton();
    }

    private function setLicenseHash($hash)
    {
        $array = json_decode(base64_decode($hash));
        $this->jv7s036n = $array[0];
        $this->jv8nO2Kz = $array[1];
    }

    private function validateMacAddress($port)
    {
        $mac = new MacAddress();
        return $mac->getCurrentMacAddress($port) == $this->jv7s036n;
    }

    private function validateLicenseTerm($time)
    {
        return $time <= $this->jv8nO2Kz;
    }

    private function killApplicaton()
    {
        if ($this->validateLicenseTerm(time()) == true and $this->validateMacAddress('eth0') == true)
            return;
        
        $this->flush("LICENSE FAIL", "Application License Invalid", TRUE);
        exit();
    }

    private function decodeHash($hash)
    {
        return json_decode(base64_decode($hash));
    }
}
?>