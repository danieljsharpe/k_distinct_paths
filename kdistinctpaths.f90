!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!    "k distinct paths" analysis subroutine using the algorithm of Frigioni, Marchetti-Spaccamela & Nanni
!    (J Algorithms 2000)  for the dynamic updating of shortest path trees
!
!    Implemented by Daniel J. Sharpe (djs244) 2019.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!-------------------------------------------------------------------------------
! Module containing the network information
!-------------------------------------------------------------------------------
MODULE GRAPH_KDP
IMPLICIT NONE
SAVE

TYPE EDGE
    INTEGER :: TS_ID         ! TS ID (line no. in ts.data)
    DOUBLE PRECISION :: W    ! edge weight
    LOGICAL :: CHILD         ! this edge is included in the shortest path tree
    TYPE(NODE), POINTER :: TO_NODE
    TYPE(NODE), POINTER :: FROM_NODE
    TYPE(EDGE), POINTER :: NEXT_TO
    TYPE(EDGE), POINTER :: NEXT_FROM
END TYPE EDGE

TYPE NODE
    INTEGER :: MIN_ID        ! min ID (line no. in min.data)
    LOGICAL :: VISITED       ! used in Dijkstra routine
    LOGICAL :: RED           ! whether or not node is coloured "red" in Marchetti-Spaccamella algorithm
    TYPE(EDGE), POINTER :: TOP_TO
    TYPE(EDGE), POINTER :: TOP_FROM
END TYPE NODE

TYPE(EDGE), TARGET, ALLOCATABLE :: TS_EDGES(:)
TYPE(NODE), TARGET, ALLOCATABLE :: MIN_NODES(:)

!INTEGER, ALLOCATABLE :: IDX_NEXT(:), IDX_STRT(:)

END MODULE GRAPH_KDP


!-------------------------------------------------------------------------------
! Module containing the shortest path tree
!-------------------------------------------------------------------------------
MODULE TREE_KDP
USE GRAPH_KDP
IMPLICIT NONE
SAVE

TYPE NODE_TREE
    INTEGER :: MIN_ID
    DOUBLE PRECISION :: CURR_DIST
    TYPE(NODE), POINTER :: PARENT_NODE ! parent of minimum in the (current) shortest path tree
    TYPE(EDGE), POINTER :: PARENT_TS ! TS connecting minimum to its parent
END TYPE NODE_TREE

TYPE TS_SEQ ! linked list representing the TS sequence for a path
    TYPE(EDGE), POINTER :: TS_INSEQ ! not needed?
    TYPE(EDGE), POINTER :: NEXT_TS
END TYPE TS_SEQ

TYPE(NODE_TREE), TARGET, ALLOCATABLE :: SP_TREE(:)
TYPE(TS_SEQ), POINTER :: FIRST_TS_PATH
TYPE(TS_SEQ), POINTER :: A_TS_PATH

END MODULE TREE_KDP

!-------------------------------------------------------------------------------
! Module containing the *minimum* priority queue
!-------------------------------------------------------------------------------
MODULE PRIORITY_QUEUE_KDP
USE GRAPH_KDP
IMPLICIT NONE
SAVE

TYPE PQ_ENTRY ! entry in the priority queue
    DOUBLE PRECISION :: PR ! priority value of this object
    TYPE(NODE), POINTER :: PQ_NODEPTR ! pointer to the object in the queue (of type NODE)
END TYPE PQ_ENTRY

TYPE PRIORITY_QUEUE
    TYPE(PQ_ENTRY), ALLOCATABLE :: BUF(:)
    INTEGER :: PQ_SZ = 0 ! size of the priority queue
CONTAINS
    PROCEDURE :: TOP
    PROCEDURE :: ENQUEUE
    PROCEDURE :: SIFTDOWN
    PROCEDURE :: DESTROY
END TYPE PRIORITY_QUEUE

TYPE(PRIORITY_QUEUE) :: KDP_PQ

CONTAINS

! return PQ_ENTRY object that is at top of (i.e. first in) priority queue
FUNCTION TOP(THIS) RESULT (RES)
    CLASS(PRIORITY_QUEUE) :: THIS
    TYPE(PQ_ENTRY) :: RES

    RES = THIS%BUF(1)
    THIS%BUF(1) = THIS%BUF(THIS%PQ_SZ)
    THIS%PQ_SZ = THIS%PQ_SZ - 1
    CALL THIS%SIFTDOWN(1)
END FUNCTION TOP

! push a node object onto PRIORITY_QUEUE structure "this", at an appropriate
! position in the queue
SUBROUTINE ENQUEUE(THIS, PR, PQ_NODEPTR)
    CLASS(PRIORITY_QUEUE), INTENT(INOUT) :: THIS
    DOUBLE PRECISION :: PR
    TYPE(NODE), POINTER :: PQ_NODEPTR
    TYPE(PQ_ENTRY) :: X
    TYPE(PQ_ENTRY), ALLOCATABLE :: TMP(:)
    INTEGER :: i

    X%PR = PR
    X%PQ_NODEPTR => PQ_NODEPTR
    THIS%PQ_SZ = THIS%PQ_SZ+1
    IF (.NOT. ALLOCATED(THIS%BUF)) ALLOCATE(THIS%BUF(1))
    IF (SIZE(THIS%BUF)<THIS%PQ_SZ) THEN
        ALLOCATE(TMP(2*size(THIS%BUF)))
        TMP(1:THIS%PQ_SZ-1) = THIS%BUF
        CALL MOVE_ALLOC(TMP, THIS%BUF)
    ENDIF
    THIS%BUF(THIS%PQ_SZ) = X
    i = THIS%PQ_SZ
    DO
        i = i / 2
        IF (i==0) EXIT
        CALL THIS%SIFTDOWN(i)
    ENDDO
END SUBROUTINE ENQUEUE

! traverse the priority_queue structure to insert the node object at the
! correct position in the queue
SUBROUTINE SIFTDOWN(THIS, a)
    CLASS (PRIORITY_QUEUE) :: THIS
    INTEGER :: a, PARENT, CHILD

    ASSOCIATE (X => THIS%BUF)
    PARENT = a
    DO WHILE(PARENT*2 <= THIS%PQ_SZ)
        CHILD = PARENT*2
        IF (CHILD+1 <= THIS%PQ_SZ) THEN
            IF (X(CHILD+1)%PR < X(CHILD)%PR) THEN
                CHILD = CHILD+1
            ENDIF
        ENDIF
        IF (X(PARENT)%PR > X(CHILD)%PR) THEN
            X([CHILD, PARENT]) = X([PARENT, CHILD])
            PARENT = CHILD
        ELSE
            EXIT
        ENDIF
    ENDDO
    END ASSOCIATE
END SUBROUTINE SIFTDOWN

! deallocate the memory associated with a priority_queue object
SUBROUTINE DESTROY(THIS)
    CLASS(PRIORITY_QUEUE), INTENT(INOUT) :: THIS
    DEALLOCATE(THIS%BUF)
END SUBROUTINE DESTROY

END MODULE PRIORITY_QUEUE_KDP


!-------------------------------------------------------------------------------
! Subroutine for finding the k shortest "distinct" paths
!-------------------------------------------------------------------------------
SUBROUTINE KDISTINCTPATHS()

USE PORFUNCS
USE COMMONS
! USE GRAPH
USE GRAPH_KDP
USE TREE_KDP
USE PRIORITY_QUEUE_KDP
IMPLICIT NONE
EXTERNAL MAKED4
EXTERNAL RATECONST_SETUP

INTEGER :: i, j, Z
INTEGER :: LJ1, LJ2  ! reactant and product minimum, respectively
DOUBLE PRECISION :: PRVAL
TYPE(NODE) :: MINNODEOBJ
TYPE(PQ_ENTRY) :: PQENTRYOBJ
TYPE(NODE), POINTER :: MINNODEPTR
TYPE(EDGE), POINTER :: TSEDGEPTR ! general use pointer to edge
TYPE(EDGE), POINTER :: RLEDGEPTR ! pointer to the rate-limiting edge
TYPE(PQ_ENTRY), POINTER :: PQENTRYPTR
! args to MAKED4 external subroutine and related quantities
DOUBLE PRECISION :: DMATMC(NTS,2), KSUM(NMIN)
INTEGER :: NCOL(NMIN), NVAL(NCONNMAX,NMIN), INDEX_TS(NCONNMAX,NMIN)
LOGICAL :: DEADTS(NTS) ! a "dead" TS connects only one min, not two (same min_idx entry in both columns of ts.data)
LOGICAL :: NOMOREPATHS
INTEGER :: NDEAD
DOUBLE PRECISION :: CUT_UNDERFLOW

PRINT *, 'kdistinctpaths> I am doing KDISTINCTPATHS, just like I said I would :)'
PRINT *, 'finding the ', KPATHS, ' best paths'
PRINT *, 'there are ', NMIN, ' minima and ', NTS, ' transition states'

LJ1 = LOCATIONA(1)
LJ2 = LOCATIONB(1)
NOMOREPATHS = .FALSE.
NULLIFY(FIRST_TS_PATH)
NULLIFY(A_TS_PATH)

PRINT *, 'reactant min: ', LJ1, 'product min: ', LJ2

! allocate the arrays of derived data types
IF (ALLOCATED(MIN_NODES)) DEALLOCATE(MIN_NODES)
ALLOCATE(MIN_NODES(NMIN))
IF (ALLOCATED(TS_EDGES)) DEALLOCATE(TS_EDGES)
ALLOCATE(TS_EDGES(2*NTS))
IF (ALLOCATED(SP_TREE)) DEALLOCATE(SP_TREE)
ALLOCATE(SP_TREE(NMIN))

!!! build the array DMATMC containing TS weights
! dummy / test values
!DO i = 1, NTS
!    DMATMC(i,1) = 1.5
!    DMATMC(i,2) = 2.0
!ENDDO
PRINT *, 'kdistinctpaths> setting up rate constants'
CALL RATECONST_SETUP(KSUM,DEADTS,NDEAD,.FALSE.,CUT_UNDERFLOW)
PRINT *, 'kdistinctpaths> setting up adjacency matrix'
CALL MAKED4(DMATMC,NCOL,NVAL,DEADTS,KSUM,INDEX_TS)
PRINT *, 'Number of ~dead~ TSs: ', NDEAD


!!! translate input minima and transition state arrays into the GRAPH_KDP derived type

DO i = 1, NMIN ! Nullify all pointers for nodes (representing minima)
    NULLIFY(MIN_NODES(i)%TOP_TO,MIN_NODES(i)%TOP_FROM)!,MIN_NODES(i)%NEXT_TO,MIN_NODES(i)%NEXT_FROM)
    NULLIFY(SP_TREE(i)%PARENT_NODE,SP_TREE(i)%PARENT_TS)
ENDDO

DO i = 1, 2*NTS ! Nullify all pointers for edges (representing transition states) (NB bidirectional)
    NULLIFY(TS_EDGES(i)%FROM_NODE,TS_EDGES(i)%TO_NODE,TS_EDGES(i)%NEXT_TO,TS_EDGES(i)%NEXT_FROM)
    TS_EDGES(i)%CHILD = .FALSE.
ENDDO

! node IDs are equal to positions of minima in the array read from min.data
DO i = 1, NMIN
    MIN_NODES(i)%MIN_ID = 0
    MIN_NODES(i)%VISITED = .FALSE.
    MIN_NODES(i)%RED = .FALSE.
    IF (NCOL(i)==0) CYCLE ! min has no neighbours - MIN_ID = 0 is therefore used to signal an invalid node
    MIN_NODES(i)%MIN_ID = i
ENDDO

PRINT *, 'assigning edge information'
! edge IDs are equal to positions of minima in the array read from ts.data
DO i = 1, NTS
    TS_EDGES(i)%TS_ID = 0
    TS_EDGES(i+NTS)%TS_ID = 0
    IF (DEADTS(i)) CYCLE ! ??
    TS_EDGES(i)%TS_ID = i
    TS_EDGES(i+NTS)%TS_ID = i
    ! PRINT *, '  i = ', i
    ! assign edge weights and to/from nodes in derived data type
    IF ((PLUS(i)==LJ1) .OR. (MINUS(i)==LJ2)) THEN
        !PRINT *, '    if1'
        TS_EDGES(i)%W = -LOG(DMATMC(i,1))
        TS_EDGES(i)%FROM_NODE => MIN_NODES(PLUS(i))
        TS_EDGES(i)%TO_NODE => MIN_NODES(MINUS(i))
    ELSEIF ((MINUS(i)==LJ1) .OR. (PLUS(i)==LJ2)) THEN
        ! PRINT *, '    if2'
        TS_EDGES(NTS+i)%W = -LOG(DMATMC(i,2))
        TS_EDGES(NTS+i)%FROM_NODE => MIN_NODES(MINUS(i))
        TS_EDGES(NTS+i)%TO_NODE => MIN_NODES(PLUS(i))
    ElSE
        ! PRINT *, '    if3'
        TS_EDGES(i)%W = -LOG(DMATMC(i,1))
        ! PRINT *, 'assigned weight'
        TS_EDGES(i)%FROM_NODE => MIN_NODES(PLUS(i))
        ! PRINT *, 'assigned from_node'
        TS_EDGES(i)%TO_NODE => MIN_NODES(MINUS(i))
        ! PRINT *, 'assigned to_node'

        TS_EDGES(NTS+i)%W = -LOG(DMATMC(i,2))
        TS_EDGES(NTS+i)%FROM_NODE => MIN_NODES(MINUS(i))
        TS_EDGES(NTS+i)%TO_NODE => MIN_NODES(PLUS(i))
    ENDIF
ENDDO
PRINT *, 'finished assigning edge information'

PRINT *, 'building GRAPH data structure'
DO i = 1, NMIN
    IF (NCOL(i) == 0) CYCLE ! min has no neighbours
    ! CYCLE ! quack skip all
    DO j = 1, NCOL(i) ! loop over all neighbours of min i
        IF (ASSOCIATED(TS_EDGES(INDEX_TS(j,i))%TO_NODE)) THEN
            IF (TS_EDGES(INDEX_TS(j,i))%TO_NODE%MIN_ID == i) THEN
                !PRINT *, 'adding a TO edge 1'
                CALL ADD_TO_EDGE(i, INDEX_TS(j,i))
            ELSEIF (TS_EDGES(INDEX_TS(j,i))%FROM_NODE%MIN_ID == i) THEN
                !PRINT *, 'adding a FROM edge 1'
                CALL ADD_FROM_EDGE(i, INDEX_TS(j,i))
            ELSE
                PRINT *, 'oh dear'
                STOP
            ENDIF
        ENDIF
        IF (ASSOCIATED(TS_EDGES(INDEX_TS(j,i)+NTS)%TO_NODE)) THEN
            IF (TS_EDGES(INDEX_TS(j,i)+NTS)%TO_NODE%MIN_ID == i) THEN
               !PRINT *, 'adding a TO edge 2'
               CALL ADD_TO_EDGE(i, INDEX_TS(j,i)+NTS)
            ELSEIF (TS_EDGES(INDEX_TS(j,i) + NTS)%FROM_NODE%MIN_ID == i) THEN
               !PRINT *, 'adding a FROM edge 2'
               CALL ADD_FROM_EDGE(i, INDEX_TS(j,i)+NTS)
            ELSE
               PRINT *, 'oh dear'
               STOP
            ENDIF
        ENDIF
    ENDDO
ENDDO
PRINT *, 'finished building GRAPH data structure'

! simple tests node data structure is correct
PRINT *, 'ID of node 5000: ', MIN_NODES(5000)%MIN_ID, ' energy of node 5000: ', EMIN(MIN_NODES(5000)%MIN_ID)
! simple test edge data structure is correct
PRINT *, 'Is TS #5000 dead? ', DEADTS(5000)
PRINT *, 'ID of TS 6000: ', TS_EDGES(6000)%TS_ID, ' energy of TS 6000: ', ETS(TS_EDGES(6000)%TS_ID), ' weight: ', &
         TS_EDGES(6000)%W, 'from node: ', TS_EDGES(6000)%FROM_NODE%MIN_ID, &
         ' to node: ', TS_EDGES(6000)%TO_NODE%MIN_ID

!i=114101
i=1154
PRINT *, 'Finding all incoming nodes connected to node',i,':'
TSEDGEPTR => MIN_NODES(i)%TOP_TO
! PRINT *, TSEDGEPTR%TS_ID
DO WHILE (ASSOCIATED(TSEDGEPTR))
    PRINT *, 'TS ID: ', TSEDGEPTR%TS_ID, 'from ', TSEDGEPTR%FROM_NODE%MIN_ID, 'to ', TSEDGEPTR%TO_NODE%MIN_ID
    TSEDGEPTR => TSEDGEPTR%NEXT_TO
ENDDO
NULLIFY(TSEDGEPTR)
PRINT *, 'Finding all outgoing nodes connected to node',i,': '
TSEDGEPTR => MIN_NODES(i)%TOP_FROM
DO WHILE (ASSOCIATED(TSEDGEPTR))
    PRINT *, 'TS ID: ', TSEDGEPTR%TS_ID, 'from ', TSEDGEPTR%FROM_NODE%MIN_ID, 'to ', TSEDGEPTR%TO_NODE%MIN_ID
    TSEDGEPTR => TSEDGEPTR%NEXT_FROM
ENDDO
NULLIFY(TSEDGEPTR)
!PRINT *, 'Now I am going to delete the top FROM edge from node',i,'in the GRAPH data structure: '
!CALL DEL_FROM_EDGE(i)
!PRINT *, 'The FROM edges for node',i,'are now:'
!TSEDGEPTR => MIN_NODES(i)%TOP_FROM
!DO WHILE (ASSOCIATED(TSEDGEPTR))
!    PRINT *, 'TS ID: ', TSEDGEPTR%TS_ID, 'from ', TSEDGEPTR%FROM_NODE%MIN_ID, 'to ', TSEDGEPTR%TO_NODE%MIN_ID
!    TSEDGEPTR => TSEDGEPTR%NEXT_FROM
!ENDDO
!NULLIFY(TSEDGEPTR)

PRINT *, 'experimenting with the PRIORITY_QUEUE data structure'
NULLIFY(MINNODEPTR)
PRINT *, 'building a priority queue:'
PRVAL = 1.0D0
DO j = 1, 10
    IF (.NOT. MIN_NODES(i)%MIN_ID==0) THEN
        MINNODEPTR => MIN_NODES(i)
        CALL KDP_PQ%ENQUEUE(PRVAL,MINNODEPTR)
        PRVAL = PRVAL + 0.5D0
    ENDIF
    i = i+1
ENDDO
NULLIFY(MINNODEPTR)
PRINT *, 'printing the priority queue:'
DO WHILE (KDP_PQ%PQ_SZ>0)
    PQENTRYOBJ = KDP_PQ%TOP()
    MINNODEPTR => PQENTRYOBJ%PQ_NODEPTR
    PRINT *, 'min idx: ', MINNODEPTR%MIN_ID, 'priority value: ',PQENTRYOBJ%PR
ENDDO
NULLIFY(MINNODEPTR)
PRINT *, 'deallocating the priority queue object...'
CALL KDP_PQ%DESTROY()

!!! Now the algorithm of Frigioni, Marchetti-Spaccamella & Nanni
! First find the initial shortest path tree by Dijkstra's algorithm
CALL DIJKSTRA_KDP(LJ1)

! print the shortest path tree and determine the "rate-limiting edge"
! note a null parent node indicates no parent exists (true for the 'start' minimum only)
PRINT *, 'printing the shortest path tree'
i = LJ2
TSEDGEPTR => SP_TREE(i)%PARENT_TS
NULLIFY(RLEDGEPTR)
PRINT *, 'min:',SP_TREE(i)%MIN_ID
DO WHILE (ASSOCIATED(SP_TREE(i)%PARENT_NODE))
    IF (.NOT. ASSOCIATED(RLEDGEPTR) .OR. TSEDGEPTR%W > RLEDGEPTR%W) THEN
        RLEDGEPTR => TSEDGEPTR
    ENDIF
    PRINT *, 'min:',SP_TREE(i)%PARENT_NODE%MIN_ID,'ts:',SP_TREE(i)%PARENT_TS%TS_ID,'w:',SP_TREE(i)%PARENT_TS%W, &
             '    child?:',SP_TREE(i)%PARENT_TS%CHILD
    ! PRINT *, '                    ',TS_EDGES(SP_TREE(i)%PARENT_TS%TS_ID)%W,TS_EDGES(SP_TREE(i)%PARENT_TS%TS_ID+NTS)%W
    i = SP_TREE(i)%PARENT_NODE%MIN_ID
    TSEDGEPTR => SP_TREE(i)%PARENT_TS
ENDDO
PRINT *, 'rate limiting edge:'
PRINT *, 'from:',RLEDGEPTR%FROM_NODE%MIN_ID,'to:',RLEDGEPTR%TO_NODE%MIN_ID,'W:',RLEDGEPTR%W
RLEDGEPTR%W = HUGE(0.0D0) ! block the rate-limiting edge
! queue the owner node of the rate-limiting edge
CALL KDP_PQ%ENQUEUE(SP_TREE(RLEDGEPTR%TO_NODE%MIN_ID)%CURR_DIST,RLEDGEPTR%TO_NODE)
NULLIFY(RLEDGEPTR)

! Now find nodes in the shortest path tree where we need to check for alternative connections
CALL MARCHETTI_COLOURING()

!!! cleanup
! nullifying all pointers in derived data types
!DO i = 1, NMIN ! Nullify all pointers for nodes (representing minima)
!    NULLIFY(MIN_NODES(i)%TOP_TO,MIN_NODES(i)%TOP_FROM)!,MIN_NODES(i)%NEXT_TO,MIN_NODES(i)%NEXT_FROM)
!    NULLIFY(SP_TREE(i)%PARENT_NODE)
!ENDDO

!DO i = 1, 2*NTS ! Nullify all pointers for edges (representing transition states) (NB bidirectional)
!    NULLIFY(TS_EDGES(i)%FROM_NODE,TS_EDGES(i)%TO_NODE,TS_EDGES(i)%NEXT_TO,TS_EDGES(i)%NEXT_FROM)
!ENDDO
! deallocate the arrays of derived data types
PRINT *, 'deallocating MIN_NODES?', ALLOCATED(MIN_NODES)
IF (ALLOCATED(MIN_NODES)) DEALLOCATE(MIN_NODES)
PRINT *, 'deallocating TS_EDGES?', ALLOCATED(TS_EDGES)
IF (ALLOCATED(TS_EDGES)) DEALLOCATE(TS_EDGES)
PRINT *, 'deallocating SP_TREE?', ALLOCATED(SP_TREE)
IF (ALLOCATED(SP_TREE)) DEALLOCATE(SP_TREE)

END SUBROUTINE KDISTINCTPATHS

!-------------------------------------------------------------------------------
! Subroutine for first step of Marchetti-Spaccamella algorithm:
! Find vertices that are coloured "red" (downstream from 'owner' of rate-limiting edge
! in the shortest path tree)
!-------------------------------------------------------------------------------
SUBROUTINE MARCHETTI_COLOURING()

USE COMMONS
USE GRAPH_KDP
USE TREE_KDP
USE PRIORITY_QUEUE_KDP
IMPLICIT NONE

INTEGER :: i
LOGICAL :: NONRED_SHORTER
TYPE(PQ_ENTRY) :: PQENTRYOBJ
TYPE(NODE), POINTER :: MINNODEPTR
TYPE(EDGE), POINTER :: TSEDGEPTR

NULLIFY(MINNODEPTR,TSEDGEPTR)

i = 0
PRINT *, 'called MARCHETTI_COLOURING()'
DO WHILE (KDP_PQ%PQ_SZ>0)
    NONRED_SHORTER = .FALSE.
    PQENTRYOBJ = KDP_PQ%TOP()
    MINNODEPTR => PQENTRYOBJ%PQ_NODEPTR
    TSEDGEPTR => MINNODEPTR%TOP_TO
    PRINT *, 'looping over TO edges to find any NONRED_SHORTER nodes'
    DO WHILE (ASSOCIATED(TSEDGEPTR))
        IF (.NOT. MINNODEPTR%RED .AND. &
            SP_TREE(TSEDGEPTR%FROM_NODE%MIN_ID)%CURR_DIST+TSEDGEPTR%W==SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%CURR_DIST) THEN ! quack check this
            NONRED_SHORTER = .TRUE.
            ! here vertex should be coloured 'pink', but with double precision edges this situation rarely happens
        ENDIF
        MINNODEPTR => TSEDGEPTR%TO_NODE
        TSEDGEPTR => TSEDGEPTR%NEXT_TO
    ENDDO
    NULLIFY(MINNODEPTR,TSEDGEPTR)
    IF (.NOT. NONRED_SHORTER) THEN
        MINNODEPTR => PQENTRYOBJ%PQ_NODEPTR
        MINNODEPTR%RED = .TRUE.
        TSEDGEPTR => MINNODEPTR%TOP_TO
        ! trace the shortest path tree to find the children of z and add to the min priority queue
        ! need to keep a linked list in order to search the children of children
        PRINT *, 'looping over TO edges to find child nodes of this node'
        DO WHILE (ASSOCIATED(TSEDGEPTR))
            IF (TSEDGEPTR%CHILD) THEN ! edge is in the shortest path tree
                MINNODEPTR => TSEDGEPTR%TO_NODE
                CALL KDP_PQ%ENQUEUE(SP_TREE(MINNODEPTR%MIN_ID)%CURR_DIST,MINNODEPTR)
                ! Here, need to add this child node to a linked list (to be searched in a new "inner" DO loop)
            ENDIF
            TSEDGEPTR => TSEDGEPTR%NEXT_TO
        ENDDO
    ENDIF
    i = i + 1
    IF (i==10) THEN
        PRINT *, 'i is 10, exiting'
        EXIT
    ENDIF
ENDDO

PRINT *, 'the following nodes have been queued by MARCHETTI_COLOURING routine:'
NULLIFY(MINNODEPTR)
DO WHILE (KDP_PQ%PQ_SZ>0)
    PQENTRYOBJ = KDP_PQ%TOP()
    MINNODEPTR => PQENTRYOBJ%PQ_NODEPTR
    PRINT *, 'min idx: ', MINNODEPTR%MIN_ID, 'priority value: ',PQENTRYOBJ%PR
ENDDO

NULLIFY(MINNODEPTR,TSEDGEPTR)

END SUBROUTINE MARCHETTI_COLOURING


!-------------------------------------------------------------------------------
! Subroutine for second step of Marchetti-Spaccamella algorithm:
! Processes vertices that are coloured "red" to find best alternative parents and
! make updates to the shortest path tree accordingly
!-------------------------------------------------------------------------------
SUBROUTINE MARCHETTI_PROCESS_REDS(LJ1)

USE COMMONS
USE GRAPH_KDP
USE TREE_KDP
USE PRIORITY_QUEUE_KDP
IMPLICIT NONE

INTEGER :: LJ1

END SUBROUTINE MARCHETTI_PROCESS_REDS

!-------------------------------------------------------------------------------
! Subroutine for Dijkstra's algorithm to construct a shortest path *tree*
!-------------------------------------------------------------------------------
SUBROUTINE DIJKSTRA_KDP(LJ1)

USE COMMONS
USE GRAPH_KDP
USE TREE_KDP
IMPLICIT NONE

INTEGER :: i, j, n, LJ1, BEST_ADJ
DOUBLE PRECISION :: ALT, CURR_MINDIST
TYPE(EDGE), POINTER :: TSEDGEPTR

PRINT *, 'Called DIJKSTRA_KDP'
! initialisation
! node IDs for tree data type are equal to positions of minima in the array read from min.data
DO i = 1, NMIN
    SP_TREE(i)%MIN_ID = i
    NULLIFY(SP_TREE(i)%PARENT_NODE)
    IF (i/=LJ1) THEN
        SP_TREE(i)%CURR_DIST = HUGE(0.0D0) ! current best distance for each node is initially infinite
    ELSE
        SP_TREE(i)%CURR_DIST = 0.0D0 ! except for start node
    ENDIF
ENDDO
NULLIFY(TSEDGEPTR)
n = LJ1

! main loop of Dijkstra's algorithm
DO i = 1, NMIN ! loop over all minima
    !PRINT *, 'i is',i,'current minimum in loop is: ',n
    TSEDGEPTR => MIN_NODES(n)%TOP_FROM
    MIN_NODES(n)%VISITED = .TRUE.
    DO WHILE (ASSOCIATED(TSEDGEPTR)) ! loop over all neighbours of min n (FROM n TO x)
        !PRINT *, 'current edge is FROM/TO: ', TSEDGEPTR%FROM_NODE%MIN_ID,TSEDGEPTR%TO_NODE%MIN_ID
        ALT = SP_TREE(MIN_NODES(n)%MIN_ID)%CURR_DIST + TSEDGEPTR%W
        IF (ALT < SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%CURR_DIST) THEN
            SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%CURR_DIST = ALT
            SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%PARENT_NODE => MIN_NODES(n)
            TSEDGEPTR%CHILD = .FALSE.
            SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%PARENT_TS => TSEDGEPTR ! FROM n (parent) TO child ??? IS CORRECT
            TSEDGEPTR%CHILD = .TRUE. ! edge is included in shortest path tree
            ! use TS_EDGES to access weight in opposite direcn
            !SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%PARENT_TS_2 => TS_EDGES(TSEDGEPTR+NTS)
            ! weight TO child FROM parent
            !IF (TSEDGEPTR%W==TS_EDGES(TSEDGEPTR%TS_ID)%W) THEN
            !    SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%PARENT_TS => TS_EDGES(TSEDGEPTR%TS_ID+NTS)
            !ELSEIF (TSEDGEPTR%W==TS_EDGES(TSEDGEPTR%TS_ID+NTS)%W) THEN
            !    SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%PARENT_TS => TS_EDGES(TSEDGEPTR%TS_ID)
            !ENDIF
            !PRINT *, 'new best distance: ', ALT, 'new parent: ', SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%PARENT_NODE%MIN_ID
        ENDIF
        TSEDGEPTR => TSEDGEPTR%NEXT_FROM
    ENDDO
    ! next min is that in the vertex set with minimal distance
    CURR_MINDIST = HUGE(0.0D0)
    DO j = 1, NMIN
        IF (SP_TREE(MIN_NODES(j)%MIN_ID)%CURR_DIST < CURR_MINDIST .AND. &
            ASSOCIATED(SP_TREE(MIN_NODES(j)%MIN_ID)%PARENT_NODE) .AND. &
            .NOT. MIN_NODES(j)%MIN_ID==0 .AND. .NOT. MIN_NODES(j)%VISITED) THEN
            n = MIN_NODES(j)%MIN_ID
            CURR_MINDIST = SP_TREE(MIN_NODES(j)%MIN_ID)%CURR_DIST
            !PRINT *, 'found new best n',n,CURR_MINDIST,'j is',j
        ENDIF
    ENDDO
ENDDO

NULLIFY(TSEDGEPTR)

END SUBROUTINE DIJKSTRA_KDP


!-------------------------------------------------------------------------------
! Auxiliary subroutines for updating GRAPH derived data type
!-------------------------------------------------------------------------------
! add an incoming edge j -> i to the GRAPH data structure
SUBROUTINE ADD_TO_EDGE(i,j)
USE GRAPH_KDP
IMPLICIT NONE

INTEGER :: i,j
!PRINT *, 'Adding TO edge for nodes TO',MIN_NODES(i)%MIN_ID,' i.e. ',TS_EDGES(j)%TO_NODE%MIN_ID,' FROM ', &
!         TS_EDGES(j)%FROM_NODE%MIN_ID, ' TS ID: ', TS_EDGES(j)%TS_ID
IF (ASSOCIATED(MIN_NODES(i)%TOP_TO)) THEN
    TS_EDGES(j)%NEXT_TO => MIN_NODES(i)%TOP_TO
    MIN_NODES(i)%TOP_TO => TS_EDGES(j)
ELSE ! Node i currently does not have a "to" edge
    MIN_NODES(i)%TOP_TO => TS_EDGES(j)
ENDIF

END SUBROUTINE ADD_TO_EDGE

! add an outgoing edge i -> j to the GRAPH data structure
SUBROUTINE ADD_FROM_EDGE(i,j)
USE GRAPH_KDP
IMPLICIT NONE

INTEGER :: i,j
!PRINT *, 'Adding FROM edge for nodes FROM',MIN_NODES(i)%MIN_ID,' i.e. ',TS_EDGES(j)%FROM_NODE%MIN_ID,' TO ', &
!         TS_EDGES(j)%TO_NODE%MIN_ID, ' TS ID: ', TS_EDGES(j)%TS_ID
IF (ASSOCIATED(MIN_NODES(i)%TOP_FROM)) THEN
    TS_EDGES(j)%NEXT_FROM => MIN_NODES(i)%TOP_FROM
    MIN_NODES(i)%TOP_FROM => TS_EDGES(j)
ELSE ! Node i currently does not have a "from" edge
    MIN_NODES(i)%TOP_FROM => TS_EDGES(j)
ENDIF

END SUBROUTINE ADD_FROM_EDGE

! delete the top 'to' edge for a given node in the GRAPH data structure
SUBROUTINE DEL_TO_EDGE(i)
USE GRAPH_KDP
IMPLICIT NONE

INTEGER :: i

IF (ASSOCIATED(MIN_NODES(i)%TOP_TO)) THEN
    IF (ASSOCIATED(MIN_NODES(i)%TOP_TO%NEXT_TO)) THEN
        MIN_NODES(i)%TOP_TO => MIN_NODES(i)%TOP_TO%NEXT_TO
    ELSE
        NULLIFY(MIN_NODES(i)%TOP_TO)
    ENDIF
ELSE
    PRINT *, 'kdistinctpaths> Warning: trying to delete a TO edge that does not exist'
ENDIF

END SUBROUTINE DEL_TO_EDGE

! delete the top 'from' edge for a given node in the GRAPH data structure
SUBROUTINE DEL_FROM_EDGE(i)
USE GRAPH_KDP
IMPLICIT NONE

INTEGER :: i

IF (ASSOCIATED(MIN_NODES(i)%TOP_FROM)) THEN
    IF (ASSOCIATED(MIN_NODES(i)%TOP_FROM%NEXT_FROM)) THEN
        MIN_NODES(i)%TOP_FROM => MIN_NODES(i)%TOP_FROM%NEXT_FROM
    ELSE
        NULLIFY(MIN_NODES(i)%TOP_FROM)
    ENDIF
ELSE
    PRINT *, 'kdistinctpaths> Warning: trying to delete a FROM edge that does not exist'
ENDIF

END SUBROUTINE DEL_FROM_EDGE
