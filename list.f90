module linked_list

  implicit none

  type::list
     !! Doubly-linked list
     type(list_node),pointer::first => null()
         !! First node of list
     type(list_node),pointer::last => null()
         !! Last node of list
   contains
     procedure::push_front
     procedure::push_back
     procedure::remove
  end type list

  type::list_node
     !! Node of a doubly-linked list
     class(*),pointer::data => null()
         !! Node data
     type(list),pointer::parent => null()
         !! Pointer to the list
     type(list_node),pointer::prev => null()
         !! Pointer to left child
     type(list_node),pointer::next => null()
         !! Pointer to right child
   contains
     procedure::insert_before
     procedure::insert_after
  end type list_node

contains

  subroutine push_front(this,data)
    !! Insert `data` at the front of list `this`

    ! Arguments
    class(list),intent(inout),target::this
        !! List object
    class(*),intent(in),target::data

    type(list_node),pointer::new_node
        !! New list node

    allocate(new_node)
    
    new_node%data=>data
    new_node%parent=>this

    if(associated(this%first)) then
       this%first%prev=>new_node
       new_node%next=>this%first
    else
       this%last=>new_node
    end if

    this%first=>new_node
       
  end subroutine push_front

  subroutine push_back(this,data)
    !! Insert `data` at the back of list `this`

    ! Arguments
    class(list),intent(inout),target::this
        !! List object
    class(*),intent(in),target::data
        !! Data to insert

    type(list_node),pointer::new_node

    allocate(new_node)

    new_node%data=>data
    new_node%parent=>this

    if(associated(this%last)) then
       this%last%next=>new_node
       new_node%prev=>this%last
    else
       this%first=>new_node
    end if

    this%last=>new_node
       
  end subroutine push_back

  subroutine insert_before(this,data)
    !! Insert `data` before list node `this`

    ! Arguments
    class(list_node),intent(inout),target::this
        !! List object
    class(*),intent(in),target::data
        !! Data to insert

    type(list_node),pointer::new_node
        !! New list node

    allocate(new_node)

    new_node%data=>data

    if(associated(this%prev)) then
       new_node%prev=>this%prev
       this%prev%next=>new_node
    else
       this%parent%first=>new_node
       new_node%prev=>null()
    end if

    new_node%next=>this
    this%prev=>new_node

  end subroutine insert_before

  subroutine insert_after(this,data)
    !! Insert `data` after list node `this`

    ! Arguments
    class(list_node),intent(inout),target::this
        !! List object
    class(*),intent(in),target::data
        !! Data to insert

    type(list_node),pointer::new_node
        !! New list node

    allocate(new_node)

    new_node%data=>data

    if(associated(this%prev)) then
       new_node%prev=>this
    else
       this%parent%last=>new_node
       new_node%next=>null()
    end if

    this%next=>new_node

  end subroutine insert_after

  subroutine remove(this,delete_node)
    !! Remove node from list

    ! Arguments
    class(list),intent(inout)::this
    type(list_node),intent(inout),pointer::delete_node
        !! List object

    if(associated(delete_node%prev)) then
       delete_node%prev%next=>delete_node%next
    else
       this%first=>delete_node%next
    end if

    if(associated(delete_node%next)) then
       delete_node%next%prev=>delete_node%prev
    else
       this%last=>delete_node%prev
    end if

  end subroutine remove

end module linked_list
