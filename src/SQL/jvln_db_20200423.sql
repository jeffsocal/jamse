CREATE DATABASE  IF NOT EXISTS `jvln` /*!40100 DEFAULT CHARACTER SET utf8 */;
USE `jvln`;
-- MySQL dump 10.13  Distrib 5.7.17, for Win64 (x86_64)
--
-- Host: localhost    Database: jvln
-- ------------------------------------------------------
-- Server version	5.7.29-0ubuntu0.18.04.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `accesses`
--

DROP TABLE IF EXISTS `accesses`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `accesses` (
  `pk` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `user_pk` int(10) unsigned NOT NULL,
  `file_pk` int(10) unsigned NOT NULL,
  `owner` tinyint(4) NOT NULL DEFAULT '0',
  `read` tinyint(4) NOT NULL DEFAULT '0',
  `write` tinyint(4) NOT NULL DEFAULT '0',
  `delete` tinyint(4) NOT NULL DEFAULT '0',
  PRIMARY KEY (`pk`),
  UNIQUE KEY `pk_UNIQUE` (`pk`),
  KEY `accesses_users_fk_idx` (`user_pk`),
  KEY `accesses_files_fk_idx` (`file_pk`),
  CONSTRAINT `accesses_files_fk` FOREIGN KEY (`file_pk`) REFERENCES `files` (`pk`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `accesses_users_fk` FOREIGN KEY (`user_pk`) REFERENCES `users` (`pk`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `credentials`
--

DROP TABLE IF EXISTS `credentials`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `credentials` (
  `pk` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `login` varchar(255) DEFAULT NULL,
  `password_hashed` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`pk`),
  UNIQUE KEY `pk_UNIQUE` (`pk`),
  UNIQUE KEY `login_UNIQUE` (`login`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `experiments`
--

DROP TABLE IF EXISTS `experiments`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `experiments` (
  `pk` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `file_pk` int(10) unsigned NOT NULL,
  `experiment_name` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`pk`),
  KEY `file_pk` (`file_pk`,`experiment_name`),
  CONSTRAINT `exp_file_pk_fk` FOREIGN KEY (`file_pk`) REFERENCES `files` (`pk`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `files`
--

DROP TABLE IF EXISTS `files`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `files` (
  `pk` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `file_name` varchar(255) NOT NULL,
  `file_id` varchar(45) DEFAULT NULL,
  `path` varchar(255) DEFAULT NULL,
  `scans_n` int(11) DEFAULT NULL,
  PRIMARY KEY (`pk`),
  UNIQUE KEY `pk_UNIQUE` (`pk`),
  KEY `file_index` (`file_name`)
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `methods`
--

DROP TABLE IF EXISTS `methods`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `methods` (
  `pk` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `parameters` longblob NOT NULL,
  `timestamp` timestamp NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`pk`),
  UNIQUE KEY `pk_UNIQUE` (`pk`)
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `metrics`
--

DROP TABLE IF EXISTS `metrics`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `metrics` (
  `pk` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `file_pk` int(10) unsigned NOT NULL,
  `metric` varchar(45) NOT NULL,
  `n` int(11) DEFAULT NULL,
  `max` float DEFAULT NULL,
  `median` float DEFAULT NULL,
  `mean` float DEFAULT NULL,
  `min` float DEFAULT NULL,
  `q25` float DEFAULT NULL,
  `q75` float DEFAULT NULL,
  `q005` float DEFAULT NULL,
  `q995` float DEFAULT NULL,
  PRIMARY KEY (`pk`),
  KEY `file_pk` (`file_pk`),
  KEY `metric` (`metric`),
  CONSTRAINT `file_pk_fk` FOREIGN KEY (`file_pk`) REFERENCES `files` (`pk`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `psm_alts`
--

DROP TABLE IF EXISTS `psm_alts`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `psm_alts` (
  `pk` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `scan_pk` int(10) unsigned NOT NULL DEFAULT '0',
  `peptide` varchar(255) NOT NULL DEFAULT 'na',
  `pre_mma` float NOT NULL DEFAULT '0',
  `psm_mma` float NOT NULL DEFAULT '0',
  `psm_z` float unsigned DEFAULT '0',
  `psm_i` int(10) unsigned DEFAULT '0',
  `psm_o` int(10) unsigned DEFAULT '0',
  `psm_score` float unsigned NOT NULL DEFAULT '0',
  `psm_fcorr` float unsigned NOT NULL DEFAULT '0',
  `psm_prob` float unsigned NOT NULL DEFAULT '0',
  `hit_prob` float unsigned NOT NULL DEFAULT '0',
  `alt_prob` float unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`pk`),
  UNIQUE KEY `pk_UNIQUE` (`pk`),
  KEY `peptide_index` (`peptide`),
  KEY `psm_score_index` (`psm_prob`,`pre_mma`),
  KEY `scan_pk_fk0` (`scan_pk`),
  CONSTRAINT `scan_pk_fk0` FOREIGN KEY (`scan_pk`) REFERENCES `scans` (`pk`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `psms`
--

DROP TABLE IF EXISTS `psms`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `psms` (
  `pk` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `scan_pk` int(10) unsigned NOT NULL DEFAULT '0',
  `method_pk` int(10) unsigned NOT NULL DEFAULT '0',
  `peptide` varchar(255) NOT NULL DEFAULT 'na',
  `pre_mma` float NOT NULL DEFAULT '0',
  `psm_mma` float NOT NULL DEFAULT '0',
  `psm_z` float unsigned DEFAULT '0',
  `psm_i` int(10) unsigned DEFAULT '0',
  `psm_o` int(10) unsigned DEFAULT '0',
  `psm_score` float unsigned NOT NULL DEFAULT '0',
  `psm_fcorr` float unsigned NOT NULL DEFAULT '0',
  `psm_prob` float unsigned NOT NULL DEFAULT '0',
  `hit_prob` float unsigned NOT NULL DEFAULT '0',
  `rank_prob` float unsigned NOT NULL DEFAULT '0',
  `hit_rank` int(10) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`pk`),
  UNIQUE KEY `pk_UNIQUE` (`pk`),
  KEY `peptide_index` (`peptide`),
  KEY `psm_score_index` (`psm_prob`,`pre_mma`,`hit_rank`),
  KEY `scan_pk_fk` (`scan_pk`),
  KEY `psms_methods_pk_fk` (`method_pk`),
  CONSTRAINT `psms_methods_pk_fk` FOREIGN KEY (`method_pk`) REFERENCES `methods` (`pk`) ON DELETE CASCADE ON UPDATE NO ACTION,
  CONSTRAINT `scan_pk_fk` FOREIGN KEY (`scan_pk`) REFERENCES `scans` (`pk`) ON DELETE CASCADE ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `scans`
--

DROP TABLE IF EXISTS `scans`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `scans` (
  `pk` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `file_pk` int(10) unsigned NOT NULL,
  `method_pk` int(10) unsigned NOT NULL,
  `scan_id` int(10) unsigned NOT NULL,
  `title` varchar(255) DEFAULT NULL,
  `rt_sec` float DEFAULT NULL,
  `pre_mz` float DEFAULT NULL,
  `pre_z` int(11) DEFAULT NULL,
  `pre_int` float DEFAULT NULL,
  `ms_it_sec` float DEFAULT NULL,
  `ms_ct` varchar(16) DEFAULT NULL,
  `ms_ce` float DEFAULT NULL,
  `peaks` longblob,
  `spectrum` longblob,
  `int_rsd` float DEFAULT NULL,
  `int_rad` float DEFAULT NULL,
  `int_iqr` float DEFAULT NULL,
  `int_ntp` float DEFAULT NULL,
  `int_prd` float DEFAULT NULL,
  `t_peaks` int(11) DEFAULT NULL,
  `n_peaks` int(11) DEFAULT NULL,
  `spec_fcorr` float DEFAULT NULL,
  `spec_prob` float DEFAULT NULL,
  PRIMARY KEY (`pk`),
  UNIQUE KEY `pk_UNIQUE` (`pk`),
  KEY `precursor_index` (`rt_sec`,`pre_mz`,`pre_z`),
  KEY `quality_index` (`spec_fcorr`,`spec_prob`,`n_peaks`,`t_peaks`),
  KEY `files_pk_fk_idx` (`file_pk`),
  KEY `methods_pk_fk_idx` (`method_pk`),
  CONSTRAINT `files_pk_fk` FOREIGN KEY (`file_pk`) REFERENCES `files` (`pk`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `scan_methods_pk_fk` FOREIGN KEY (`method_pk`) REFERENCES `methods` (`pk`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Temporary view structure for view `summary_experiments`
--

DROP TABLE IF EXISTS `summary_experiments`;
/*!50001 DROP VIEW IF EXISTS `summary_experiments`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE VIEW `summary_experiments` AS SELECT 
 1 AS `experiment_name`,
 1 AS `n_files`,
 1 AS `n_scans`*/;
SET character_set_client = @saved_cs_client;

--
-- Temporary view structure for view `summary_ids`
--

DROP TABLE IF EXISTS `summary_ids`;
/*!50001 DROP VIEW IF EXISTS `summary_ids`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE VIEW `summary_ids` AS SELECT 
 1 AS `experiment_name`,
 1 AS `file_id`,
 1 AS `file_name`,
 1 AS `scans_n`,
 1 AS `scans_num`,
 1 AS `psms_num`,
 1 AS `ratio_id`*/;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `users`
--

DROP TABLE IF EXISTS `users`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `users` (
  `pk` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `cred_pk` int(10) unsigned NOT NULL,
  `name_last` varchar(255) DEFAULT NULL,
  `name_first` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`pk`),
  UNIQUE KEY `pk_UNIQUE` (`pk`),
  UNIQUE KEY `cred_pk_UNIQUE` (`cred_pk`),
  CONSTRAINT `users_cred_pk_fk` FOREIGN KEY (`cred_pk`) REFERENCES `credentials` (`pk`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `validations`
--

DROP TABLE IF EXISTS `validations`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `validations` (
  `pk` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `psm_pk` int(10) unsigned NOT NULL,
  `fdr_metric` varchar(45) NOT NULL,
  `method_pk` int(10) unsigned NOT NULL,
  `fdr_eprob` float unsigned NOT NULL,
  `fdr_prob` float unsigned DEFAULT NULL,
  `fdr_qvalue` float unsigned DEFAULT NULL,
  PRIMARY KEY (`pk`),
  UNIQUE KEY `pk_UNIQUE` (`pk`),
  KEY `perf_index` (`fdr_eprob`,`fdr_prob`,`fdr_qvalue`),
  KEY `psms_pk_fk_idx` (`psm_pk`),
  KEY `methods_pk_fk_idx` (`method_pk`),
  CONSTRAINT `pref_methods_pk_fk` FOREIGN KEY (`method_pk`) REFERENCES `methods` (`pk`) ON DELETE CASCADE ON UPDATE NO ACTION,
  CONSTRAINT `psms_pk_fk` FOREIGN KEY (`psm_pk`) REFERENCES `psms` (`pk`) ON DELETE CASCADE ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Final view structure for view `summary_experiments`
--

/*!50001 DROP VIEW IF EXISTS `summary_experiments`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = cp850 */;
/*!50001 SET character_set_results     = cp850 */;
/*!50001 SET collation_connection      = cp850_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`scbi`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `summary_experiments` AS select `experiments`.`experiment_name` AS `experiment_name`,count(0) AS `n_files`,sum(`files`.`scans_n`) AS `n_scans` from (`experiments` join `files`) where (`experiments`.`file_pk` = `files`.`pk`) group by `experiments`.`experiment_name` */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `summary_ids`
--

/*!50001 DROP VIEW IF EXISTS `summary_ids`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`scbi`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `summary_ids` AS select `ms`.`experiment_name` AS `experiment_name`,`ms`.`file_id` AS `file_id`,`ms`.`file_name` AS `file_name`,`ms`.`scans_n` AS `scans_n`,`ms`.`scans_num` AS `scans_num`,`id`.`psms_num` AS `psms_num`,round((`id`.`psms_num` / `ms`.`scans_num`),3) AS `ratio_id` from (((select `jvln`.`experiments`.`experiment_name` AS `experiment_name`,`jvln`.`files`.`file_id` AS `file_id`,`jvln`.`files`.`file_name` AS `file_name`,`jvln`.`files`.`scans_n` AS `scans_n`,count(0) AS `scans_num` from ((`jvln`.`scans` join `jvln`.`files`) join `jvln`.`experiments`) where ((`jvln`.`files`.`pk` = `jvln`.`experiments`.`file_pk`) and (`jvln`.`files`.`pk` = `jvln`.`scans`.`file_pk`)) group by `jvln`.`files`.`pk`,`jvln`.`experiments`.`experiment_name`,`jvln`.`files`.`file_id`,`jvln`.`files`.`file_name`,`jvln`.`files`.`scans_n`)) `ms` left join (select `jvln`.`files`.`file_id` AS `file_id`,count(0) AS `psms_num` from ((`jvln`.`files` join `jvln`.`scans`) join `jvln`.`psms`) where ((`jvln`.`files`.`pk` = `jvln`.`scans`.`file_pk`) and (`jvln`.`scans`.`pk` = `jvln`.`psms`.`scan_pk`) and (`jvln`.`psms`.`hit_rank` = 1) and (`jvln`.`psms`.`psm_prob` >= 0.95)) group by `jvln`.`files`.`file_id`,`jvln`.`files`.`pk`) `id` on((`id`.`file_id` = `ms`.`file_id`))) */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2020-04-23 15:44:24
