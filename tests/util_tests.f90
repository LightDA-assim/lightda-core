module util_tests
  use util
  use tree
  implicit none
contains
  subroutine test_append_array
    real(kind=8),allocatable::a(:),b(:)
    real(kind=8)::x
    integer,parameter::n=10
    integer::i

    allocate(a(n))

    do i=1,n
       a(i)=i
    end do

    x=n+1

    b=append_array(a,x)

    do i=1,n
       if(b(i)/=a(i)) then
          print *,'Error: b(',i,')=',b(i),', expected b(',i,')=',a(i)
          error stop
       end if
    end do

    if(b(n+1)/=x) then
       print *,'Error: b(',n+1,')=',b(n+1),', expected b(',n+1,')=',x
       error stop
    end if

  end subroutine test_append_array

  subroutine check_list(test_list,data,read_inds)
    use linked_list

    type(list),intent(in)::test_list
    integer,intent(in),target::data(:)
    integer,intent(in)::read_inds(:)
    type(list_node),pointer::last_node
    type(list_node),pointer::current_node
    integer::i

    if(size(read_inds)==0 .and. associated(test_list%first)) then
       print *,'Error: First is not null but zero size list expected'
       error stop
    end if

    if(size(read_inds)==0 .and. associated(test_list%last)) then
       print *,'Error: Last is not null but zero size list expected'
       error stop
    end if

    if(size(read_inds)>0 .and. .not.associated(test_list%first)) then
       print *,'Error: First is null but nonzero list size expected'
       error stop
    end if

    if(size(read_inds)>0 .and. .not.associated(test_list%last)) then
       print *,'Error: Last is null but nonzero list size expected'
       error stop
    end if

    if(associated(test_list%first).and.associated(test_list%first%prev)) then
       print *,'Error: first%prev is not null'
       error stop
    end if

    if(associated(test_list%last).and.associated(test_list%last%next)) then
       print *,'Error: last%next is not null'
       error stop
    end if

    ! Iterate over list and check contents
    current_node=>test_list%first
    i=1
    do while(associated(current_node))
       if(i>size(read_inds)) then
          print *,'Error: Reached end of read_inds but not end of list.'
          error stop
       end if
       if(.not.associated(current_node%data,data(read_inds(i)))) then
          print *,'Error: Wrong data in node',i
          error stop
       end if
       if(i>1 .and..not.associated(current_node%prev%data,data(read_inds(i-1)))) then
          print *,'Error: Wrong prev for node ',i
          error stop
       end if
       if(i>1) then
          if(.not.associated(current_node%prev,last_node)) then
             print *,'Error: Wrong prev for node ',i
             error stop
          end if
       end if
       last_node=>current_node
       current_node=>current_node%next
       i=i+1
    end do

    if(associated(test_list%last).and..not.associated(last_node,test_list%last)) then
       print *,'Error: Wrong last node'
       error stop
    end if

  end subroutine check_list

  subroutine test_linked_list
    use linked_list

    type(list)::test_list
    type(list_node),pointer::current_node
    integer,target::data(4)
    integer::i

    ! Populate data array
    data=(/1,2,3,4/)

    call check_list(test_list,data,(/ integer :: /))

    ! Push data into list
    do i=1,4
       call test_list%push_back(data(i))
    end do

    ! Check contents of list
    call check_list(test_list,data,(/1,2,3,4/))

    ! Remove last node
    call test_list%remove(test_list%last)

    ! Check contents of list
    call check_list(test_list,data,(/1,2,3/))

    ! Push last node back on list
    call test_list%push_back(data(4))

    ! Find third node
    current_node=>test_list%first
    do i=1,2
       current_node=>current_node%next
    end do

    ! Delete third node
    call test_list%remove(current_node)

    ! Check contents of list
    call check_list(test_list,data,(/1,2,4/))

    ! Push data(3) to front of list
    call test_list%push_front(data(3))

    ! Check contents of list
    call check_list(test_list,data,(/3,1,2,4/))

    ! Delete last node
    call test_list%remove(test_list%last)

    ! Check contents of list
    call check_list(test_list,data,(/3,1,2/))

    ! Insert before last node
    call test_list%last%insert_before(data(4))

    ! Check contents of list
    call check_list(test_list,data,(/3,1,4,2/))

    ! Empty list
    do i=1,4
       call test_list%remove(test_list%last)
    end do

  end subroutine test_linked_list

  subroutine test_tree

    type(node),target::root
        !! Root node
    type(node),pointer::newRoot
        !! New root after deleting a node
    type(node),target::insertNode1,insertNode2,insertNode4,insertNode5, &
         insertNode6,insertNode7,insertNode8
        !! Nodes to be inserted into tree
    type(node),pointer::foundNode
        !! Node found in tree
    real(kind=8),target::data(8)
        !! Array to be stored in tree

    ! Fill data array
    data=(/1,2,3,4,5,6,7,8/)

    ! Assign data elements and keys to nodes

    root%key=3
    root%data=>data(3)

    insertNode1%key=1
    insertNode1%data=>data(1)

    insertNode2%key=2
    insertNode1%data=>data(2)

    insertNode4%key=4
    insertNode4%data=>data(4)

    insertNode5%key=5
    insertNode5%data=>data(5)

    insertNode6%key=6
    insertNode6%data=>data(6)

    insertNode7%key=7
    insertNode7%data=>data(7)

    insertNode8%key=8
    insertNode8%data=>data(8)

    ! Insert nodes into tree
    call root%insert(insertNode5)
    call root%insert(insertNode4)
    call root%insert(insertNode1)
    call root%insert(insertNode2)
    call root%insert(insertNode7)
    call root%insert(insertNode6)
    call root%insert(insertNode8)

    ! Test that all the nodes can be found
    foundNode=>root%find(1)
    if(.not.associated(foundNode,insertNode1)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    foundNode=>root%find(2)
    if(.not.associated(foundNode,insertNode2)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    foundNode=>root%find(3)
    if(.not.associated(foundNode,root)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    foundNode=>root%find(4)
    if(.not.associated(foundNode,insertNode4)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    foundNode=>root%find(5)
    if(.not.associated(foundNode,insertNode5)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    foundNode=>root%find(6)
    if(.not.associated(foundNode,insertNode6)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    foundNode=>root%find(7)
    if(.not.associated(foundNode,insertNode7)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    foundNode=>root%find(8)
    if(.not.associated(foundNode,insertNode8)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    ! Test that we get an unassociated pointer when searching for a key that
    ! is not in the tree
    foundNode=>root%find(-1)
    if(associated(foundNode)) then
       print *,'Error: Node returned from node%find when searching for ', &
            ' a nonexistent key. Should have returned an unassociated', &
            ' pointer instead.'
       error stop
    end if

    ! Now test deleting a node from the tree
    newRoot=>root%delete(5)

    ! Make sure that node%delete returned the correct pointer
    if(.not.associated(newRoot,root)) then
       print *,'Error: Wrong pointer returned from node%delete'
       error stop
    end if

    ! Make sure the deleted key is gone
    foundNode=>root%find(5)
    if(associated(foundNode)) then
       print *,'Error: node%delete failed to remove the specified key'
       error stop
    end if

    ! Make sure the rest of the nodes are still there
    foundNode=>root%find(1)
    if(.not.associated(foundNode,insertNode1)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    foundNode=>root%find(2)
    if(.not.associated(foundNode,insertNode2)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    foundNode=>root%find(3)
    if(.not.associated(foundNode,root)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    foundNode=>root%find(4)
    if(.not.associated(foundNode,insertNode4)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    foundNode=>root%find(6)
    if(.not.associated(foundNode,insertNode6)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    foundNode=>root%find(7)
    if(.not.associated(foundNode,insertNode7)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if

    foundNode=>root%find(8)
    if(.not.associated(foundNode,insertNode8)) then
       print *,'Error: Wrong node returned from node%find'
       error stop
    end if


  end subroutine test_tree
end module util_tests

program test_util
  use util_tests
  implicit none

  call test_append_array
  call test_linked_list
  call test_tree
end program test_util
