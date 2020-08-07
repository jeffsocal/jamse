CREATE PROCEDURE `psms_by_fileid` (OUT param_01 char(32))
BEGIN
SELECT psms.pk pk, 
	  scan_id, rt_sec, pre_mz, pre_z, 
	  int_rsd, int_rad, int_iqr, int_ntp, 
	  spec_fcorr, spec_prob, 
	  pre_mma, psm_z, psm_i, psm_o, psm_mma, 
	  psm_score, psm_fcorr, 
	  psm_prob, hit_prob, rank_prob, hit_rank
	  FROM jvln.psms, jvln.scans, jvln.files
	  WHERE jvln.scans.file_pk = jvln.files.pk 
	  AND jvln.psms.scan_pk = jvln.scans.pk 
	  AND file_id = param_01;
END
