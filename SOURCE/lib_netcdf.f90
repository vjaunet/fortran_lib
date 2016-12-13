MODULE lib_netcdf
  use netcdf
  implicit none
  !
  !
  !
  !        netCDF interface for handling data
  !
  ! Description:
  !         wrapper to  netcdf fortran 90 routines to handle all the hard work
  !         simple to the aboc routines will open, read and store data
  !         in dedicated objects.
  !
  ! author :
  !        vincent jaunet
  !
  ! licence : GPL v3.0
  !
  !============================================================================
  ! data containers
  type ncdf_var_t
     character(len=20)                       :: name
     real(kind=8), dimension(:), allocatable :: data
  end type ncdf_var_t
  type ncdf_dim_t
     character(len=20)                       :: name
     integer                                 :: len
  end type ncdf_dim_t

  type netcdf_data
     !netcdf variables------------------------------------
     ! file id
     integer                                 ::ncdf_id=-1

     ! dimension handlers
     integer                                 ::ndim,nvar,nattrib
     integer                                 ::unlimDimId, formatnum

     !data containers
     type(ncdf_dim_t), dimension(:), allocatable ::dimensions
     type(ncdf_var_t),  dimension(:), allocatable ::var

     ! type bounded procedures
   contains
     procedure          :: create => lib_netcdf_create
     procedure          :: destroy => lib_netcdf_destroy
     procedure          :: getvarid => lib_netcdf_getvarid
     procedure          :: setVarNames => lib_netcdf_setVarNames
     procedure          :: setDimNames => lib_netcdf_setDimNames
     procedure          :: fillVar => lib_netcdf_fillVar
     procedure          :: fillDim => lib_netcdf_fillDim
     procedure          :: open_file => lib_netcdf_open
     procedure          :: close_file => lib_netcdf_close
     procedure          :: read_data => lib_netcdf_read
     procedure          :: write_data => lib_netcdf_write
     procedure          :: print_help => lib_netcdf_print_help
     procedure          :: print_info

  end type netcdf_data


contains
  subroutine lib_netcdf_open(this,filename)
    class(netcdf_data)          ::this
    character(len=*)            ::filename
    !--------------------------------------

    call check(nf90_open(trim(filename),NF90_NOWRITE,this%ncdf_id))


  end subroutine lib_netcdf_open

  subroutine lib_netcdf_close(this)
    class(netcdf_data)          ::this
    !---------------------------------

    call check(nf90_close(this%ncdf_id))

    this%ncdf_id = -1

  end subroutine lib_netcdf_close

  !********************************************************
  !
  !   Constructor/Destructor Section
  !
  !
  !********************************************************
  subroutine lib_netcdf_create(this,&
       dim_len,dim_name,&
       var_len,var_name,&
       nattrib,unlimDimId,formatnum)
    class(netcdf_data)                          ::this
    integer          ,dimension(:)              ::dim_len
    character(len=*) ,dimension(:), optional    ::dim_name
    integer                                     ::var_len
    character(len=*),dimension(var_len),optional::var_name
    integer   ,optional         ::unlimDimId, formatnum, nattrib

    integer                                     ::ndim,nvar
    integer                                     ::total_size

    integer                                     ::ivar,idim
    !-------------------------------------------------------------

    this%ndim=size(dim_len)
    this%nvar=var_len

    this%unlimDimId=0
    this%nattrib=0
    this%formatnum=3

    if (present(unlimDimId)) this%unlimDimId=unlimDimId
    if (present(nattrib))    this%nattrib=nattrib
    if (present(formatnum))  this%formatnum=formatnum

    !dim container memory allocation
    allocate(this%dimensions(this%ndim))
    do idim=1,this%ndim
       this%dimensions(idim)%len=dim_len(idim)
    end do

    !-- data container memory allocation
    allocate(this%var(this%nvar))
    !for the dimensions
    do ivar=1,this%ndim
       allocate(this%var(ivar)%data(this%dimensions(idim)%len))
    end do
    !for the variables
    total_size=product(dim_len(1:this%ndim))
    do ivar=this%ndim+1,this%nvar
       allocate(this%var(ivar)%data(total_size))
    end do

    !set data names
    if (present(dim_name)) then
       do idim=1,this%ndim
          this%dimensions(idim)%name=trim(dim_name(idim))
       end do
    else
       do idim=1,this%ndim
          write(this%dimensions(idim)%name,'(a1,i1)')'x',idim
       end do
    end if

    if (present(var_name)) then
       do ivar=1,this%nvar
          this%var(ivar)%name=trim(var_name(ivar))
       end do
    else
       do ivar=1,this%nvar
          write(this%var(ivar)%name,'(a1,i1)')'v',ivar
       end do
    end if

  end subroutine lib_netcdf_create

  subroutine lib_netcdf_fillDim(this,dim_id,x)
    class(netcdf_data)                          ::this
    real(kind=8)       ,dimension(:)            ::x
    integer                                     ::dim_id
    !------------------------------------------------------
    if (.not.allocated(this%var)) then
       write(06,*)'lib_netcdf_fillDim : the current object has not been allocated yet'
       STOP 2
    end if

    if (dim_id > this%ndim) then
       write(06,*)'lib_netcdf_fillDim : dimension too high'
       STOP 2
    end if

    if (size(x) /= (this%dimensions(dim_id)%len)) then
       write(06,*)'lib_netcdf_fillDim : data not of the correct size'
       STOP 2
    end if


    this%var(dim_id)%data=x

  end subroutine lib_netcdf_fillDim

  subroutine lib_netcdf_fillVar(this,var_id,x)
    class(netcdf_data)                          ::this
    real(kind=8)       ,dimension(:)            ::x
    integer                                     ::var_id
    !------------------------------------------------------

    if (.not.allocated(this%var)) then
       write(06,*)'lib_netcdf_fillVar : the current object has not been allocated yet'
       STOP 2
    end if

    if (var_id <= this%ndim .or. var_id > this%nvar) then
       write(06,*)'lib_netcdf_fillVar : var_id not on the variable range'
       STOP 2
    end if

    if (size(x) /= product(this%dimensions(:)%len)) then
       write(06,*)'lib_netcdf_fillVar : data not of the correct size'
       STOP 2
    end if

    this%var(var_id)%data=x

  end subroutine lib_netcdf_fillVar

  subroutine lib_netcdf_setVarNames(this,var_name)
    class(netcdf_data)                          ::this
    character(len=*) ,dimension(:)              ::var_name
    integer                                     ::ivar
    !------------------------------------------------------
    if (.not.allocated(this%var)) then
       write(06,*)'lib_netcdf_setVarNames : the current object has not been allocated yet'
       STOP 2
    end if

    do ivar=1,this%nvar
       write(this%var(ivar)%name,'(a1,i1)')'v',ivar
    end do
  end subroutine lib_netcdf_setVarNames

  subroutine lib_netcdf_setDimNames(this,dim_name)
    class(netcdf_data)                          ::this
    character(len=*) ,dimension(:)              ::dim_name
    integer                                     ::idim
    !------------------------------------------------------

    if (.not.allocated(this%dimensions)) then
       write(06,*)'lib_netcdf_setDimNames : the current object has not been allocated yet'
       STOP 2
    end if

    do idim=1,this%ndim
       write(this%dimensions(idim)%name,'(a1,i1)')'x',idim
    end do

  end subroutine lib_netcdf_setDimNames

  subroutine lib_netcdf_destroy(this)
    class(netcdf_data)                          ::this
    integer                                     ::ivar
    !------------------------------------------------------

    if (allocated(this%dimensions)) deallocate(this%dimensions)
    do ivar=1,this%nvar
       if (allocated(this%var(ivar)%data)) deallocate(this%var(ivar)%data)
    end do
    if (allocated(this%var)) deallocate(this%var)

  end subroutine lib_netcdf_destroy


  !********************************************************
  !
  !> function getVariableid
  !> return the index of the desired variable
  !> @param charater(len=*) varname
  !> @return integer*4 varid
  !>
  !********************************************************
  integer function lib_netcdf_getvarid(this,varname)
    class(netcdf_data)                          ::this
    character(len=*)                            ::varname
    integer                                     ::ivar
    !------------------------------------------------------

    do ivar=1,this%nvar
       if (this%var(ivar)%name == varname) then
          lib_netcdf_getvarid = ivar
          return
       end if
    end do

    !the name has not been found :
    write(06,*)'lib_netcdf_getvarid : unkwnown variable name "',varname,'"'


  end function  lib_netcdf_getvarid



  !********************************************************
  !
  !   Read data section
  !
  !
  !********************************************************
  subroutine lib_netcdf_read(this,filename)
    class(netcdf_data)                 ::this
    character(len=*)   ,optional       ::filename
    integer                            ::var_id, iDim
    integer                            ::var_size=1
    integer, dimension(:), allocatable ::start_, count_,dim_id
    integer                            ::ndim
    !------------------------------------------------

    if (this%ncdf_id == -1) then
       if (present(filename)) then
          call this%open_file(filename)
       else
          STOP 'lib_netcdf_open : no filename given'
       end if
    end if

    !get to know what is in the file
    call check(nf90_inquire(this%ncdf_id,nDimensions=this%ndim,&
         nVariables=this%nvar,&
         nAttributes=this%nattrib,&
         unlimitedDimId=this%unlimDimId,&
         formatNum=this%formatnum))

    allocate(this%var(this%nvar))
    allocate(this%dimensions(this%ndim))
    do iDim=1,this%nDim
       call check(nf90_inquire_dimension(this%ncdf_id,iDim,&
            name=this%dimensions(iDim)%name,&
            len=this%dimensions(idim)%len))
    end do

    !for all the var get their name and values
    !not that we read the whole file in one step
    allocate(dim_id(this%ndim))
    dim_id=0
    do var_id=1,this%nvar

       call check(nf90_inquire_variable(this%ncdf_id,var_id,&
            name=this%var(var_id)%name, dimids=dim_id,&
            ndims=ndim))

       var_size=product(this%dimensions(dim_id(1:ndim))%len)
       allocate(this%var(var_id)%data(var_size))

       allocate(start_(ndim),count_(ndim))
       start_=1
       count_=this%dimensions(dim_id(1:ndim))%len
       call check(nf90_get_var(this%ncdf_id,&
            varid=var_id,&
            values=this%var(var_id)%data, start=start_, count=count_))

       deallocate(start_, count_)
    end do

    !close the file
    call this%close_file()

  end subroutine lib_netcdf_read

  !********************************************************
  !
  !   Write data section
  !
  !
  !********************************************************

  subroutine lib_netcdf_write(this,filename)
    class(netcdf_data)               ::this
    character(len=*)                 ::filename

    integer  ,dimension(this%ndim)   ::dim_id,dim_len
    integer  ,dimension(this%nvar)   ::var_id
    integer                          ::idim,ivar,trash
    integer                          ::var_size
    integer  ,dimension(this%ndim)   ::start_,count_
    !-------------------------------------------

    !create the data file
    call check(nf90_create(filename, NF90_CLOBBER, this%ncdf_id))

    !definitions------------------

    !Define the dimensions.
    do idim=1,this%ndim
       call check(nf90_def_dim(this%ncdf_id, &
            this%dimensions(idim)%name,&
            this%dimensions(idim)%len,&
            dim_id(idim)))
    end do

    !Define coordinate data variables
    do idim=1,this%ndim
       call check(nf90_def_var(this%ncdf_id,&
            this%var(idim)%name, NF90_DOUBLE, &
            dim_id(idim),var_id(idim)))
    end do

    !Define the other data variable
    !beware that we are in column major mode
    if (this%ndim >= 2) dim_id(1:2)=dim_id(2:1:-1)
    do ivar=this%ndim+1,this%nvar
       call check(nf90_def_var(this%ncdf_id,&
            this%var(ivar)%name, NF90_DOUBLE, &
            dim_id, var_id(ivar)))
    end do

    call check(nf90_enddef(this%ncdf_id))
    !End Definitions----------------

    !Write Data---------------------
    do ivar=1,this%nvar

       call check(nf90_inquire_variable(this%ncdf_id,&
            varid=var_id(ivar),ndims=idim,dimids=dim_id))

       start_=1
       count_=this%dimensions(dim_id(1:this%ndim))%len
       call check(nf90_put_var(this%ncdf_id,&
            varid=var_id(ivar),&
            values=this%var(ivar)%data,&
            start=start_, count=count_))
    end do

    !close the file-----------------
    call check(nf90_close(this%ncdf_id))

  end subroutine lib_netcdf_write


  subroutine print_info(this)
    class(netcdf_data)             ::this
    integer                        ::ivar,idim
    !------------------------------------------
    write(06,*)'netCDF file info :'
    write(06,*)'======================'
    write(06,*)
    write(06,*)'Dimensions :'
    write(06,*)'------------'
    do idim=1,this%ndim
       write(06,'(a,i3,a,a,a,i7)')' Dimension ',idim,' name :',trim(this%dimensions(idim)%name),&
            &', length = ', this%dimensions(idim)%len
    end do
    write(06,*)
    write(06,*)'Variables :'
    write(06,*)'------------'
    do ivar=1,this%nvar
       write(06,'(a,i02,a,a)')' Var ',ivar,' name : ',trim(this%var(ivar)%name)
    end do
    write(06,*)
  end subroutine print_info

  subroutine lib_netcdf_print_help(this)
    class(netcdf_data)             ::this
    !------------------------------------------
    write(06,*)'lib_netCDF routines :'
    write(06,*)'======================'
    write(06,*)''
    write(06,*)'-- obj.open(filename)'
    write(06,*)'-- obj.close()'
    write(06,*)'-- obj.read() : reads and store the data into '
    write(06,*)'obj.dimensions and obj.var derived type'
    write(06,*)'-- obj.print_info(), obj.print_help()'
    STOP
  end subroutine lib_netcdf_print_help


  ! Private routine retruning the exit status
  ! of called netCDF function
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
       write(06,*)trim(nf90_strerror(status))
       stop 2
    end if
  end subroutine check

END MODULE lib_netcdf
