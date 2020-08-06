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

          if(this%left%key > new_node%key) then

             ! Insert into left subtree
             call this%left%insert(new_node)

          else

             ! Insert new_node between `this` and left subtree
             new_node%right=>this%left
             this%left=>new_node

          end if

       else

          ! Insert new_node as left node of `this`
          this%left=>new_node

       end if

    else

       if(associated(this%right)) then

          if(this%right%key < new_node%key) then

             ! Insert into right subtree
             call this%right%insert(new_node)

          else

             ! Insert new_node between `this` and right subtree
             new_node%left=>this%right
             this%right=>new_node

          end if

       else

          ! Insert new_node as right node of `this`
          this%right=>new_node

       end if

    end if

  end subroutine insert

end module tree
