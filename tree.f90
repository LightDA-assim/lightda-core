module tree

  implicit none

  type :: node
     !! Node of a binary tree
     integer,pointer::key
         !! Node key
     class(*),pointer::data
         !! Node data
     type(node),pointer::left
         !! Pointer to left child
     type(node),pointer::right
         !! Pointer to right child
   contains
     procedure::find
     procedure::insert
     procedure::delete
  end type node

contains

  recursive function find(this,key) result(foundNode)
    !! Search node `this` and its children for a node whose key is equal
    !! to `key`. Returns a pointer to the node whose key is equal to `key`,
    !! or an unassociated pointer if no such node exists.

    ! Arguments
    class(node),intent(in),target::this
        !! Node to search
    integer,intent(in)::key
        !! Key to search for
    type(node),pointer::foundNode
        !! Node matching key, or an unassociated pointer if no match was found

    if(this%key==key) then

       ! We're there! Point to ourselves and exit
       foundNode=>this
       return

    elseif(key<this%key) then

       if(associated(this%left)) then

          ! Search left subtree and return
          foundNode=>this%left%find(key)
          return

       end if

    else

       if(associated(this%right)) then

          ! Search right subtree and return
          foundNode=>this%right%find(key)
          return

       end if

    end if

    ! Nothing was found return an unassociated pointer
    nullify(foundNode)

  end function find

  recursive subroutine insert(this,new_node)
    !! Insert node `new_node` into the tree below `this`.
    !!
    !! If the key of `new_node` is equal to that of an existing node,
    !! do nothing.

    ! Arguments
    class(node),intent(inout),target::this
        !! Root node of subtree
    type(node),intent(inout),target::new_node
        !! New node to insert

    if(this%key==new_node%key) then

       ! Node already exists
       return

    elseif(new_node%key<this%key) then

       if(associated(this%left)) then

          ! Insert into left subtree
          call this%left%insert(new_node)

       else

          ! Insert new_node as left node of `this`
          this%left=>new_node

       end if

    else

       if(associated(this%right)) then

          ! Insert into right subtree
          call this%right%insert(new_node)

       else

          ! Insert new_node as right node of `this`
          this%right=>new_node

       end if

    end if

  end subroutine insert

  recursive function delete(this,key) result(newRoot)
    !! Delete the node whose key is equal to `key`, if any such node exists,
    !! from the subtree below the node `this`.
    !!
    !! Returns a pointer to thenew root of the subtree. Return value will point
    !! to `this` if `this` matches `key`, otherwise it will be a pointer to
    !! one of the two immediate children of `this`.

    ! Arguments
    class(node),intent(inout),target::this
        !! Root node of subtree
    integer::key
        !! Key of node to delete

    type(node),pointer::newRoot
        !! New root of subtree

    if(key==this%key) then
       if(associated(this%left)) then
          newRoot=>this%left
          call newRoot%insert(this%right)
       else
          newRoot=>this%right
          call newRoot%insert(this%left)
       end if
    else
       newRoot=>this
       if(key>this%key) then
          this%right=>this%right%delete(key)
       else
          this%left=>this%left%delete(key)
       end if
    end if

  end function delete

end module tree
