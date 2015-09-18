; $Id$

common colors,r_orig,g_orig,b_orig,r_curr,g_curr,b_curr

common my_colors,       $
	                clr_blk_idx, $
			clr_wht_idx, $
			color_order, $
                        clr_tbl

common cbar,            $
                        cbar_chr_sz, $
                        cbar_fnt, $
                        cbar_fmt, $
                        cbar_idx, $
                        cbar_lbl_sz, $
                        cbar_lgn, $
                        cbar_psn, $
                        cbar_txt_clr, $
                        cbar_unit, $
                        cbar_clr_nbr
