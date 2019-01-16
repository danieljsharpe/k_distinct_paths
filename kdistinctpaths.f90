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
    TYPE(NODE), POINTER :: TO_NODE
    TYPE(NODE), POINTER :: FROM_NODE
    TYPE(EDGE), POINTER :: NEXT_TO
    TYPE(EDGE), POINTER :: NEXT_FROM
END TYPE EDGE

TYPE NODE
    INTEGER :: MIN_ID        ! min ID (line no. in min.data)
    LOGICAL :: PERMANENT     ! used in Dijkstra routine
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
IMPLICIT NONE
SAVE

TYPE NODE_TREE
    INTEGER :: MIN_ID
    DOUBLE PRECISION :: CURR_DIST
    TYPE(NODE_TREE), POINTER :: PARENT_NODE ! parent of minimum in the (current) shortest path tree
END TYPE NODE_TREE

TYPE(NODE_TREE), TARGET, ALLOCATABLE :: SP_TREE(:)

END MODULE TREE_KDP


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

INTEGER :: i, j
INTEGER :: LJ1, LJ2  ! reactant and product minimum, respectively
TYPE(EDGE), POINTER :: TSEDGEPTR
! args to MAKED4 external subroutine and related quantities
DOUBLE PRECISION :: DMATMC(NTS,2), KSUM(NMIN)
INTEGER :: NCOL(NMIN), NVAL(NCONNMAX,NMIN), INDEX_TS(NCONNMAX,NMIN)
LOGICAL :: DEADTS(NTS) ! a "dead" TS connects only one min, not two (same min_idx entry in both columns of ts.data)
INTEGER :: NDEAD
DOUBLE PRECISION :: CUT_UNDERFLOW

PRINT *, 'kdistinctpaths> I am doing KDISTINCTPATHS, just like I said I would :)'
PRINT *, 'finding the ', KPATHS, ' best paths'
PRINT *, 'there are ', NMIN, ' minima and ', NTS, ' transition states'

LJ1 = LOCATIONA(1)
LJ2 = LOCATIONB(1)

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
    NULLIFY(SP_TREE(i)%PARENT_NODE)
ENDDO

DO i = 1, 2*NTS ! Nullify all pointers for edges (representing transition states) (NB bidirectional)
    NULLIFY(TS_EDGES(i)%FROM_NODE,TS_EDGES(i)%TO_NODE,TS_EDGES(i)%NEXT_TO,TS_EDGES(i)%NEXT_FROM)
ENDDO

! node IDs are equal to positions of minima in the array read from min.data
DO i = 1, NMIN
    MIN_NODES(i)%MIN_ID = 0
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
                PRINT *, 'adding a TO edge 1'
                CALL ADD_TO_EDGE(i, INDEX_TS(j,i))
            ELSEIF (TS_EDGES(INDEX_TS(j,i))%FROM_NODE%MIN_ID == i) THEN
                PRINT *, 'adding a FROM edge 1'
                CALL ADD_FROM_EDGE(i, INDEX_TS(j,i))
            ELSE
                PRINT *, 'oh dear'
                STOP
            ENDIF
        ENDIF
        IF (ASSOCIATED(TS_EDGES(INDEX_TS(j,i)+NTS)%TO_NODE)) THEN
            IF (TS_EDGES(INDEX_TS(j,i)+NTS)%TO_NODE%MIN_ID == i) THEN
               PRINT *, 'adding a TO edge 2'
               CALL ADD_TO_EDGE(i, INDEX_TS(j,i)+NTS)
            ELSEIF (TS_EDGES(INDEX_TS(j,i) + NTS)%FROM_NODE%MIN_ID == i) THEN
               PRINT *, 'adding a FROM edge 2'
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

i=114101
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
PRINT *, 'Now I am going to delete the top TO edge from node',i,'in the GRAPH data structure: '
!CALL DEL_TO_EDGE(i)
PRINT *, 'The TO edges for node',i,'are now:'
TSEDGEPTR => MIN_NODES(i)%TOP_TO
DO WHILE (ASSOCIATED(TSEDGEPTR))
    PRINT *, 'TS ID: ', TSEDGEPTR%TS_ID, 'from ', TSEDGEPTR%FROM_NODE%MIN_ID, 'to ', TSEDGEPTR%TO_NODE%MIN_ID
    TSEDGEPTR => TSEDGEPTR%NEXT_TO
ENDDO
NULLIFY(TSEDGEPTR)

!!! Now the algorithm of Frigioni, Marchetti-Spaccamella & Nanni
! First find the initial shortest path tree by Dijkstra's algorithm
CALL DIJKSTRA_KDP(LJ1)

PRINT *, 'ID of node 5000 according to tree: ', SP_TREE(5000)%MIN_ID

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
! Subroutine for Dijkstra's algorithm to construct a shortest path *tree*
!-------------------------------------------------------------------------------
SUBROUTINE DIJKSTRA_KDP(LJ1)

USE COMMONS
USE TREE_KDP
IMPLICIT NONE

INTEGER :: i, LJ1

! node IDs for tree data type are equal to positions of minima in the array read from min.data
DO i = 1, NMIN
    SP_TREE(i)%MIN_ID = i
    IF (i/=LJ1) THEN
        SP_TREE(i)%CURR_DIST = HUGE(0.0D0) ! current best distance for each node is initially infinite
    ELSE
        SP_TREE(i)%CURR_DIST = 0.0D0 ! except for start node
    ENDIF
ENDDO



END SUBROUTINE DIJKSTRA_KDP


!-------------------------------------------------------------------------------
! Auxiliary subroutines for updating derived data types
!-------------------------------------------------------------------------------
! add an incoming edge j -> i to the GRAPH data structure
SUBROUTINE ADD_TO_EDGE(i,j)
USE GRAPH_KDP
IMPLICIT NONE

INTEGER :: i,j
PRINT *, 'Adding TO edge for nodes TO',MIN_NODES(i)%MIN_ID,' i.e. ',TS_EDGES(j)%TO_NODE%MIN_ID,' FROM ', &
         TS_EDGES(j)%FROM_NODE%MIN_ID, ' TS ID: ', TS_EDGES(j)%TS_ID
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
PRINT *, 'Adding FROM edge for nodes FROM',MIN_NODES(i)%MIN_ID,' i.e. ',TS_EDGES(j)%FROM_NODE%MIN_ID,' TO ', &
         TS_EDGES(j)%TO_NODE%MIN_ID, ' TS ID: ', TS_EDGES(j)%TS_ID
IF (ASSOCIATED(MIN_NODES(i)%TOP_FROM)) THEN
    TS_EDGES(j)%NEXT_FROM => MIN_NODES(i)%TOP_FROM
    MIN_NODES(i)%TOP_FROM => TS_EDGES(j)
ELSE ! Node i currently does not have a "from" edge
    MIN_NODES(i)%TOP_FROM => TS_EDGES(j)
ENDIF

END SUBROUTINE ADD_FROM_EDGE

! delete the top 'to' edge for a given node
!SUBROUTINE DEL_TO_EDGE(i)
!USE GRAPH_KDP
!IMPLICIT NONE

!INTEGER :: i

!IF (ASSOCIATED(MIN_NODES(i)%TOP_TO)) THEN
!    IF (ASSOCIATED(MIN_NODES(i)%TOP_TO%NEXT_TO)) THEN
!        MIN_NODES(i)%TOP_TO => MIN_NODES(i)%TOP_TO%NEXT_TO
!    ELSE
!        NULLIFY(MIN_NODES(i)%TOP_TO)
!    ENDIF
!ELSE
!    PRINT *, 'kdistinctpaths> Warning: trying to delete a TO edge that does not exist'
!ENDIF

!END SUBROUTINE DEL_TO_EDGE

! delete the top 'from' edge for a given node
!SUBROUTINE DEL_FROM_EDGE

!END SUBROUTINE DEL_FROM_EDGE
