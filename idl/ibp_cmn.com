; $Id$

@ibp_clr.com

common ibp_output,	$
                        chroma_idx, $
                        chroma_lst, $
                        chroma_nbr, $
                        dest_idx, $
                        dest_lst, $
                        dest_nbr, $
                        fmt_idx, $
                        fmt_lst, $
                        fmt_nbr, $
                        pll_idx, $
                        pll_lst, $
                        pll_nbr, $
                        shape_idx, $
                        shape_lst, $
                        shape_nbr
common ibp_sct,	$
                        dim_lst, $
                        dim_nbr, $
                        fgr_lst, $
                        fgr_nbr, $
                        fld_lst, $
                        fld_nbr, $
                        lev_dim_nm, $
                        pre_fld_lst, $
                        pre_fld_nbr, $
                        rgn_lst, $
                        rgn_nbr
common ibp_sngjd,   $
			auto_scl_jd, $
			background_jd, $
			cntr_lvl_nbr_jd, $
			cntr_lvl_min_jd, $
			cntr_ntv_jd, $
			cntr_ovl_jd, $
			info_jd, $
			misc_jd, $
			ryttl_jd, $
			ttl_jd, $
			unit_jd, $
			usr_cntr_jd, $
			xttl_jd, $
			yttl_jd

common ibp_datajd,     $	
			max_rng_x_jd, $
			rng_y_max_jd, $
			rng_x_min_jd, $
			rng_y_min_jd

common ibp_sng,	$
                        abb_sng, $
                        english_sng, $
                        fld_sng, $
                        hostname_sng, $
                        lat_sng, $
                        lev_sng, $
                        lon_sng, $
                        month_sng, $
                        rgn_sng, $
                        slice_sng, $
                        src_sng, $
                        sub_ttl_sng, $
                        sym_sng, $
                        ttl_sng, $
                        unit_pfx_sng, $
                        unit_sng, $
                        usr_sng, $
                        year_sng

common ibp_fgr,	$
                        bin_max, $
                        bin_min, $
                        bin_nbr, $
                        bin_sz, $
                        hi_lo_plt, $
                        hst_max, $
                        hst_min, $
                        lat_ctr, $
                        lon_ctr, $
                        map_lmt_dgr, $
                        map_set_pro, $
                        plt_rgn_nrm, $
                        ptr_psn, $
                        x_map_cbar_mrg, $
                        x_map_mrg, $
                        y_map_cbar_mrg, $
                        y_map_mrg
         
common ibp_data,        $
                        avg_data, $
                        cntr_ntt, $
                        cntr_fll_idx, $
                        cntr_lvl, $
                        cntr_lvl_lbl, $
                        cntr_lvl_nbr, $
                        cntr_ln_sty, $
                        cntr_lvl_max, $
                        cntr_lvl_min, $
                        cntr_ntv, $
                        cntr_thk, $
                        cntr_which_lbl, $
                        data, $
                        data_max, $
                        data_min, $
                        data_nbr, $
                        fld_idx, $
                        fld_nm, $
                        inverse_data, $
                        lat, $
                        lat_max, $
                        lat_max_idx, $
                        lat_min, $
                        lat_min_idx, $
                        lat_nbr, $
                        lat_wgt, $
                        lev, $
                        lev_max, $
                        lev_max_idx, $
                        lev_min, $
                        lev_min_idx, $
                        lev_nbr, $
                        lon, $
                        lon_max, $
                        lon_max_idx, $
                        lon_min, $
                        lon_min_idx, $
                        lon_nbr, $
                        max_rng_x, $
                        pre_fld, $
                        pre_fmt, $
                        rgn_idx, $
                        scl, $
                        time, $
                        time_max, $
                        time_max_idx, $
                        time_min, $
                        time_min_idx, $
                        time_nbr, $
                        x_dim, $
                        rng_y_max, $
                        y_dim

common ibp_constants,	$
                        False, $
                        True

common ibp_dbg,	$
                        dbg_txt_sz, $
                        dbg, $
                        dbg_vrb

common ibp_fls,	$
                        fgr_idx, $
                        fl_in, $
                        fl_out, $
                        src_nm, $
                        tbl_fl, $
                        time_lbl, $
                        time_slc, $
                        usr_dfn_clr_tbl, $
                        vld_avg_thr_pct, $
                        vrt_slc, $
                        vrt_slc_lbl, $
                        very_big_nbr

common ibp_msk,	$
                        msk_mode, $
                        msk_data, $
                        msk_fld_nm

common ibp_sv,	$
                        sv_abb_sng, $
                        sv_avg_data, $
                        sv_data, $
                        sv_fld_nm, $
                        sv_fld_sng, $
                        sv_lat, $
                        sv_lon, $
                        sv_data_max, $
                        sv_data_min, $
                        sv_data_nbr, $
                        sv_lat_nbr, $
                        sv_lon_nbr, $
                        sv_slc_sng, $
                        sv_src_sng, $
                        sv_sub_ttl_sng, $
                        sv_sym_sng, $
                        sv_unit_sng    

common ibp_wdg_1,     $
			Wcstm, $
                        Wavg_data, $
                        Wbreak, $
                        Wchroma, $
                        Wclr_tbl, $
                        Wdata_max, $
                        Wdata_min, $
                        Wdest, $
                        Wdifference, $
                        Wdim_base, $
                        Wdraw, $
                        Wfirst_lft_clm, $
                        Wfirst_rgt_clm, $
                        Wfl_in, $
                        Wfld, $
                        Wformat, $
                        Whst_base, $
                        Whst_gph, $
                        Wlat, $
                        Wlat_lev_base, $
                        Wlat_lev_gph, $
                        Wlat_list, $
                        Wlat_list_base, $
                        Wlat_max, $
                        Wlat_min, $
                        Wleft_clm, $
                        Wlev_max, $
                        Wlev_min, $
                        Wlon, $
                        Wlon_lat_base, $
                        Wlon_lat_gph, $
                        Wlon_lev_base, $
                        Wlon_lev_gph, $
                        Wlon_list, $
                        Wlon_list_base, $
                        Wlon_max, $
                        Wlon_min, $
                        Wmask, $
                        Wscl

common ibp_wdg_2,     $
                        Wfl_out, $
                        Woutput, $
                        Wpalette, $
                        Wquit, $
                        Wread, $
                        Wrgn, $
                        Wscat_base, $
                        Wscat_gph, $
                        Wsecond_lft_row, $
                        Wsecond_rgt_clm, $
                        Wshape, $
                        Wsv, $
                        Wtbl_fl, $
                        Wthird_lft_row, $
                        Wtime_lat_base, $
                        Wtime_lat_gph, $
                        Wtime_slc, $
                        Wtime_slc_base, $
                        Wtop_lvl, $
                        Wusr_hook, $
                        Wvalue, $
                        Wvrt_slc, $
                        Wvrt_slc_base, $
                        Wx_dim, $
                        Wy_dim, $
                        Wznl_avg_base, $
                        Wznl_avg_gph

common ibp_cstm,   $
			Wbin_max, $
			Wbin_min, $
			Wbin_sz, $
			Wcstm_base, $
			Wcstm_dfl, $
			Wcstm_done, $
			Wcstm_draw, $
			Wmax_rng_x, $
			Wmisc, $
			Wttl, $
			Wunit, $
			Wusr_cntr_lvl, $
			Wrng_x_min, $
			Wx_ttl, $
			Wrng_y_max, $
			Wrng_y_min, $
			Wy_ttl, $
			Wy_ttl_rgt, $
                        Wauto_scl, $
                        Wcntr_lvl, $
                        Wcntr_lvl_base, $
                        Wcntr_lvl_min, $
                        Wcntr_lvl_nbr, $
                        Wcntr_ntv, $
                        Wcntr_ovl, $
                        Wcntr_usr 

common ibp_anim,     $
                        Wanim, $
                        Wanim_base, $
                        Wanim_done, $
                        Wanim_lev, $
                        Wanim_load, $
                        Wanim_nbr, $
                        Wanim_play, $
                        Wanim_time, $
                        Wprojector, $
                        Wprojector_base, $
                        anim_idx, $
                        anim_nbr
common ibp_wdw,     $
                        hst_wdw, $
                        hst_wdw_x_sz, $
                        hst_wdw_y_sz, $
                        lat_lev_wdw, $
                        lat_lev_wdw_x_sz, $
                        lat_lev_wdw_y_sz, $
                        lon_lat_wdw, $
                        lon_lat_wdw_x_sz, $
                        lon_lat_wdw_y_sz, $
                        lon_lev_wdw, $
                        lon_lev_wdw_x_sz, $
                        lon_lev_wdw_y_sz, $
                        scat_wdw, $
                        scat_wdw_x_sz, $
                        scat_wdw_y_sz, $
                        time_lat_wdw, $
                        time_lat_wdw_x_sz, $
                        time_lat_wdw_y_sz, $
                        znl_avg_wdw, $
                        znl_avg_wdw_x_sz, $
                        znl_avg_wdw_y_sz
