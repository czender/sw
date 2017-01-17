! $Id$

! Purpose: Library of interpolation/extrapolation routines

! Copyright (C) 1994--2017 Charlie Zender
! You may copy, distribute, and/or modify this software under the terms of the GNU General Public License (GPL) Version 2
! See http://www.gnu.org/copyleft/gpl.html for full license text

! Usage: 
!use xtr_mdl ! [mdl] Extrapolation/interpolation handling

module xtr_mdl ! [mdl] Extrapolation/interpolation handling
  implicit none
  public::xtr_ini ! [sbr] Set-up extrapolation/interpolation handling

  ! Extrapolation types are as follows:
  integer,parameter,private::xtr_msk_nbr=11 ! [nbr] Number of extrapolation/interpolation masks

  logical::xtr_flg_LHS(0:xtr_msk_nbr-1) ! [enm] Extrapolation/interpolation flags
  logical::xtr_flg_RHS(0:xtr_msk_nbr-1) ! [enm] Extrapolation/interpolation flags
  character(80)::xtr_sng_LHS ! [sng] Extrapolation/interpolation string
  character(80)::xtr_sng_RHS ! [sng] Extrapolation/interpolation string

  ! nil ngh lnr are mutually exclusive
  ! frc implies prt 
  ! wgt implies frc
  
  ! NB: Keys must be unique
  integer,parameter::xtr_prt_msk=0 ! Partial extrapolation is allowed 
  integer,parameter::xtr_fll_msk=1 ! Full extrapolation is allowed (implies xtr_prt_msk)
  integer,parameter::xtr_fll_nil_msk=2 ! Set extrapolated value to 0.0
  integer,parameter::xtr_fll_ngh_msk=3 ! Set extrapolated value to value of nearest valid neighbor
  integer,parameter::xtr_fll_lnr_msk=4 ! Perform linear extrapolation using two nearest valid neighbors
  integer,parameter::xtr_prt_nil_msk=5 ! Set extrapolated value to 0.0
  integer,parameter::xtr_prt_ngh_msk=6 ! Set extrapolated value to value of nearest valid neighbor
  integer,parameter::xtr_prt_lnr_msk=7 ! Perform linear extrapolation using two nearest valid neighbors
  integer,parameter::xtr_prt_frc_msk=8 ! Use average value of overlap region
  integer,parameter::xtr_prt_wgt_msk=9 ! Set extrapolated value to average value of overlap region weighted by size of overlap region plus 0.0 weighted by size of non overlap region (implies xtr_frc_msk)
  integer,parameter::xtr_vrb_msk=10 ! Print verbose warnings whenever extrapolation is performed
  
  ! Metaflags
  integer,parameter::xtr_err=0 ! No extrapolation allowed under any circumstances
  integer,parameter::xtr_vrb=2**xtr_vrb_msk ! Print circumstances of all extrapolations

  ! Full extrapolation types. Selecting one of these sets xtr_fll_msk=.true.
  integer,parameter::xtr_fll_nil=2**xtr_fll_nil_msk ! Set extrapolated value to 0.0
  integer,parameter::xtr_fll_ngh=2**xtr_fll_ngh_msk ! Set extrapolated value to value of nearest valid neighbor
  integer,parameter::xtr_fll_lnr=2**xtr_fll_lnr_msk ! Linearly extrapolate values using values two nearest valid neighbors

  ! Partial extrapolation types. Selecting one of these sets xtr_prt_msk=.true.
  integer,parameter::xtr_prt_nil=2**xtr_prt_nil_msk ! Set extrapolated value to 0.0
  integer,parameter::xtr_prt_lnr=2**xtr_prt_lnr_msk ! Linearly extrapolate value
  integer,parameter::xtr_prt_ngh=2**xtr_prt_ngh_msk ! Set extrapolated value to value of nearest valid neighbor
  integer,parameter::xtr_prt_frc=2**xtr_prt_frc_msk ! Compute average in overlap region. Use xtr_prt_wgt_msk to determine how to weight final value
  integer,parameter::xtr_prt_wgt=2**xtr_prt_wgt_msk ! Set extrapolated value to fraction in valid region weighted by size of valid region plus 0.0 weighted by size of non-overlap region

contains

  subroutine xtr_ini( & ! [sbr] Set-up extrapolation/interpolation handling
       xtr_typ_LHS, & ! I [enm] Extrapolation type
       xtr_typ_RHS) ! I [enm] Extrapolation type
    ! Purpose: Set extrapolation flags and strings for given extrapolation types
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen,ftn_strcpy,ftn_strini,ftn_strcat ! [mdl] String manipulation
    implicit none
    ! Commons
    ! Input
    integer,intent(in)::xtr_typ_LHS ! I [enm] Extrapolation type
    integer,intent(in)::xtr_typ_RHS ! I [enm] Extrapolation type
    ! Local
    character(9) xtr_msk_sng(0:xtr_msk_nbr-1) 
    integer idx ! [idx] Counting index
    integer itmp
    integer rmn ! Remainder
    integer mlt ! Multiple (greatest factor s.t. mlt*msk<=rmn)
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Entering xtr_ini()'
    
    if (xtr_typ_LHS /= xtr_typ_RHS) then
       write (6,'(a,i4,a,i4)') &
            'xtr_ini(): WARNING xtr_typ_LHS = ',xtr_typ_LHS,'  /= xtr_typ_RHS = ',xtr_typ_RHS
    endif
    
    ! Initialize strings
    do idx=0,xtr_msk_nbr-1 ! NB: idx starts with 0
       call ftn_strini(xtr_msk_sng(idx))
       if (dbg_lvl == dbg_old) then
          write (6,'(3(a,i4),2a)') "idx = ",idx, &
               ", len(xtr_msk_sng(idx)) = ",len(xtr_msk_sng(idx)), &
               ", xtr_msk_nbr = ",xtr_msk_nbr, &
               ", xtr_msk_sng(idx) = ",xtr_msk_sng(idx)
       end if ! end if dbg
    end do ! end loop over sng
    call ftn_strcpy(xtr_msk_sng(xtr_prt_msk),'prt ')
    call ftn_strcpy(xtr_msk_sng(xtr_fll_msk),'fll ')
    call ftn_strcpy(xtr_msk_sng(xtr_fll_nil_msk),'fll_nil ')
    call ftn_strcpy(xtr_msk_sng(xtr_fll_ngh_msk),'fll_ngh ')
    call ftn_strcpy(xtr_msk_sng(xtr_fll_lnr_msk),'fll_lnr ')
    call ftn_strcpy(xtr_msk_sng(xtr_prt_nil_msk),'prt_nil ')
    call ftn_strcpy(xtr_msk_sng(xtr_prt_ngh_msk),'prt_ngh ')
    call ftn_strcpy(xtr_msk_sng(xtr_prt_lnr_msk),'prt_lnr ')
    call ftn_strcpy(xtr_msk_sng(xtr_prt_frc_msk),'prt_frc ')
    call ftn_strcpy(xtr_msk_sng(xtr_prt_wgt_msk),'prt_wgt ')
    call ftn_strcpy(xtr_msk_sng(xtr_vrb_msk),'vrb ')
    
    call ftn_strini(xtr_sng_LHS)
    call ftn_strini(xtr_sng_RHS)
    
    ! Proceed in reverse order by mask size
    rmn=xtr_typ_LHS ! Initialize remainder to original flag
    do idx=xtr_msk_nbr-1,0,-1 ! NB: idx ends with 0
       xtr_flg_LHS(idx)=.false.
       itmp=2**idx
       mlt=rmn/itmp
       rmn=mod(rmn,itmp) ! Next iteration remainder is what is leftover
       if (mlt == 1) xtr_flg_LHS(idx)=.true.
    end do ! end loop over msk
    
    ! Implement and check logical consequences of flags
    if (xtr_flg_LHS(xtr_fll_nil_msk).or.xtr_flg_LHS(xtr_fll_ngh_msk).or.xtr_flg_LHS(xtr_fll_lnr_msk)) then
       xtr_flg_LHS(xtr_fll_msk)=.true.
    else
       xtr_flg_LHS(xtr_fll_msk)=.false.
    endif ! endif
    if (xtr_flg_LHS(xtr_prt_nil_msk).or.xtr_flg_LHS(xtr_prt_ngh_msk).or. &
         xtr_flg_LHS(xtr_prt_lnr_msk).or.xtr_flg_LHS(xtr_prt_frc_msk).or. &
         xtr_flg_LHS(xtr_prt_wgt_msk)) then
       xtr_flg_LHS(xtr_prt_msk)=.true.
    else
       xtr_flg_LHS(xtr_prt_msk)=.false.
    endif ! endif
    if (xtr_flg_LHS(xtr_prt_wgt_msk)) xtr_flg_LHS(xtr_prt_frc_msk)=.true.
    
    ! Assign strings
    do idx=0,xtr_msk_nbr-1 ! NB: idx starts with 0
       if (xtr_flg_LHS(idx)) call ftn_strcat(xtr_sng_LHS,xtr_msk_sng(idx)) 
    end do ! end loop over msk
    
    ! Proceed in reverse order by mask size
    rmn=xtr_typ_RHS ! Initialize remainder to original flag
    do idx=xtr_msk_nbr-1,0,-1 ! NB: idx ends with 0
       xtr_flg_RHS(idx)=.false.
       itmp=2**idx
       mlt=rmn/itmp
       rmn=mod(rmn,itmp) ! Next iteration remainder is what is leftover
       if (mlt == 1) xtr_flg_RHS(idx)=.true.
    end do ! end loop over msk
    
    ! Implement and check logical consequences of flags
    if (xtr_flg_RHS(xtr_fll_nil_msk).or.xtr_flg_RHS(xtr_fll_ngh_msk).or.xtr_flg_RHS(xtr_fll_lnr_msk)) then
       xtr_flg_RHS(xtr_fll_msk)=.true.
    else
       xtr_flg_RHS(xtr_fll_msk)=.false.
    endif ! endif
    if (xtr_flg_RHS(xtr_prt_nil_msk).or.xtr_flg_RHS(xtr_prt_ngh_msk).or. &
         xtr_flg_RHS(xtr_prt_lnr_msk).or.xtr_flg_RHS(xtr_prt_frc_msk).or. &
         xtr_flg_RHS(xtr_prt_wgt_msk)) then
       xtr_flg_RHS(xtr_prt_msk)=.true.
    else
       xtr_flg_RHS(xtr_prt_msk)=.false.
    endif ! endif
    if (xtr_flg_RHS(xtr_prt_wgt_msk)) xtr_flg_RHS(xtr_prt_frc_msk)=.true.
    
    ! Assign strings
    do idx=0,xtr_msk_nbr-1 ! NB: idx starts with 0
       if (xtr_flg_RHS(idx)) call ftn_strcat(xtr_sng_RHS,xtr_msk_sng(idx)) 
    end do ! end loop over msk
    
    if (dbg_lvl >= dbg_io) then
       write (6,'(a,i4,a,a)') 'xtr_typ_LHS = ',xtr_typ_LHS,' = ',xtr_sng_LHS(1:ftn_strlen(xtr_sng_LHS))
       write (6,'(a,i4,a,a)') 'xtr_typ_RHS = ',xtr_typ_RHS,' = ',xtr_sng_RHS(1:ftn_strlen(xtr_sng_RHS))
    endif ! endif dbg
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Exiting xtr_ini()'
    
    return
  end subroutine xtr_ini
  
end module xtr_mdl
