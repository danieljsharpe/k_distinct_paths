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
    LOGICAL :: INFASTESTALL  ! TS is found in a fastest path, is to be printed to ts.data.fastest.all
    LOGICAL :: INFASTESTCURR ! TS is found in current fastest path, is to be printed to ts.data.fastest.x
    TYPE(NODE), POINTER :: TO_NODE
    TYPE(NODE), POINTER :: FROM_NODE
    TYPE(EDGE), POINTER :: NEXT_TO
    TYPE(EDGE), POINTER :: NEXT_FROM
END TYPE EDGE

TYPE NODE
    INTEGER :: MIN_ID        ! min ID (line no. in min.data)
    LOGICAL :: VISITED       ! used in Dijkstra routine & to keep track of entries in priority queue in Marchetti algorithm
    LOGICAL :: RED           ! whether or not node is coloured "red" in Marchetti-Spaccamella algorithm
    LOGICAL :: INFASTESTALL  ! min is found in a fastest path, is to be printed to min.data.fastest.all
    LOGICAL :: INFASTESTCURR ! min is found in current fastest path, is to be printed to ts.data.fastest.x
    TYPE(EDGE), POINTER :: TOP_TO
    TYPE(EDGE), POINTER :: TOP_FROM
END TYPE NODE

TYPE(EDGE), TARGET, ALLOCATABLE :: TS_EDGES(:)
TYPE(NODE), TARGET, ALLOCATABLE :: MIN_NODES(:)

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

TYPE NODE_LIST ! linked list of nodes
    TYPE(NODE), POINTER :: NODE_INLIST
    TYPE(NODE_LIST), POINTER :: NEXTLLENTRY
END TYPE NODE_LIST

TYPE(NODE_TREE), TARGET, ALLOCATABLE :: SP_TREE(:)
TYPE(TS_SEQ), POINTER :: FIRST_TS_PATH
TYPE(TS_SEQ), POINTER :: A_TS_PATH
TYPE(NODE_LIST), POINTER :: HEADLLRED, DYNLLRED ! second ptr is used to dynamically change target during runtime

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
! Module containing subroutines for which an explicit interface is required
!-------------------------------------------------------------------------------
MODULE KDP_SUBS
USE COMMONS
USE GRAPH_KDP
USE TREE_KDP
IMPLICIT NONE
CONTAINS

!-------------------------------------------------------------------------------
! subroutine to print a path and determine the rate limiting edge
!-------------------------------------------------------------------------------
SUBROUTINE PRINT_PATH_KDP(RLEDGEPTR,LJ2,k)

INTEGER :: i, j, k, LJ2
DOUBLE PRECISION :: SUMW
TYPE(EDGE), POINTER :: TSEDGEPTR, RLEDGEPTR
CHARACTER(LEN=8) :: FMTDESC ! format descriptor
CHARACTER(LEN=4) :: X1 ! intermediate for converting int to string using an 'internal file'
CHARACTER(LEN=11) :: FNAME
LOGICAL :: FILE_EXIST

IF (ASSOCIATED(RLEDGEPTR)) NULLIFY(RLEDGEPTR)
NULLIFY(TSEDGEPTR)
FMTDESC = '(I4.4)'
WRITE (X1,FMTDESC) k
FNAME = "Epath."//TRIM(X1)
OPEN(6,FILE=FNAME)

SUMW = 0.0D0
i = LJ2 ! note a null parent node indicates no parent exists (true for the 'start' minimum only)
j = 2
TSEDGEPTR => SP_TREE(i)%PARENT_TS
PRINT *, ''
PRINT *, 'PRINTING PATH NO.', k
PRINT *, 'min:',SP_TREE(i)%MIN_ID
WRITE(6,*) j-1, EMIN(i), SP_TREE(i)%MIN_ID
MIN_NODES(LJ2)%INFASTESTCURR = .TRUE.
DO WHILE (ASSOCIATED(SP_TREE(i)%PARENT_NODE))
    TSEDGEPTR%INFASTESTALL = .TRUE.
    SP_TREE(i)%PARENT_NODE%INFASTESTALL = .TRUE.
    TSEDGEPTR%INFASTESTCURR = .TRUE.
    SP_TREE(i)%PARENT_NODE%INFASTESTCURR = .TRUE.
    SUMW = SUMW + TSEDGEPTR%W
    IF (.NOT. ASSOCIATED(RLEDGEPTR) .OR. (TSEDGEPTR%W > RLEDGEPTR%W .AND. RLEDGE_DEFN==2)) THEN
        RLEDGEPTR => TSEDGEPTR
    ELSE IF ((ETS(TSEDGEPTR%TS_ID)-EMIN(TSEDGEPTR%FROM_NODE%MIN_ID) > &
              ETS(RLEDGEPTR%TS_ID)-EMIN(RLEDGEPTR%FROM_NODE%MIN_ID)) &
              .AND. RLEDGE_DEFN==1) THEN
        RLEDGEPTR => TSEDGEPTR ! RLEDGEPTR has max barrier height
    ELSE IF ((ETS(TSEDGEPTR%TS_ID) > ETS(RLEDGEPTR%TS_ID)) .AND. RLEDGE_DEFN==0) THEN
        RLEDGEPTR => TSEDGEPTR
    ENDIF
    PRINT *, 'min:',SP_TREE(i)%PARENT_NODE%MIN_ID,'ts:',SP_TREE(i)%PARENT_TS%TS_ID,'w:',SP_TREE(i)%PARENT_TS%W
    WRITE(6,*) j, ETS(SP_TREE(i)%PARENT_TS%TS_ID), SP_TREE(i)%PARENT_TS%TS_ID
    WRITE(6,*) j+1, EMIN(SP_TREE(i)%PARENT_NODE%MIN_ID), SP_TREE(i)%PARENT_NODE%MIN_ID
    i = SP_TREE(i)%PARENT_NODE%MIN_ID
    j = j+2
    TSEDGEPTR => SP_TREE(i)%PARENT_TS
ENDDO
PRINT *, 'rate limiting edge:'
PRINT *, 'from:',RLEDGEPTR%FROM_NODE%MIN_ID,'to:',RLEDGEPTR%TO_NODE%MIN_ID,'w:',RLEDGEPTR%W
NULLIFY(TSEDGEPTR)
CLOSE(6)
PRINT *, 'sum of weights:', SUMW, 'exp of sum of weights:', EXP(SUMW)

INQUIRE(FILE="rate_lim_cut.dat",EXIST=FILE_EXIST)
IF (FILE_EXIST) THEN
    OPEN(9,FILE="rate_lim_cut.dat",STATUS="old",POSITION="append",ACTION="write")
ELSE
    OPEN(9,FILE="rate_lim_cut.dat",STATUS="new",ACTION="write")
ENDIF
WRITE(9,*) RLEDGEPTR%TS_ID, RLEDGEPTR%W, SUMW, k
CLOSE(9)

END SUBROUTINE PRINT_PATH_KDP

!-------------------------------------------------------------------------------
! subroutine to print min.data.fastest.x and ts.data.fastest.x files for current path
!-------------------------------------------------------------------------------
SUBROUTINE PRINT_INFASTEST(LJ2,k,DOALLIN)

INTEGER :: i, k, LJ2
LOGICAL, OPTIONAL :: DOALLIN
LOGICAL :: DOALL
CHARACTER(LEN=8) :: FMTDESC ! format descriptor
CHARACTER(LEN=4) :: X1 ! intermediate for converting int to string using an 'internal file'
CHARACTER(LEN=30) :: FNAME_MDF, FNAME_TDF

IF (PRESENT(DOALLIN)) THEN
    DOALL = DOALLIN
ELSE
    DOALL = .FALSE.
ENDIF

IF (.NOT. DOALL) THEN ! print stationary point indices that are along the current path
    FMTDESC = '(I4.4)'
    WRITE (X1,FMTDESC) k
    FNAME_MDF = "min.data.fastest."//TRIM(X1)
    FNAME_TDF = "ts.data.fastest."//TRIM(X1)
    OPEN(7,FILE=FNAME_MDF)
    OPEN(8,FILE=FNAME_TDF)
ELSE
    FNAME_MDF = "min.data.fastest.all"
    FNAME_TDF = "ts.data.fastest.all"
ENDIF

OPEN(7,FILE=FNAME_MDF)
OPEN(8,FILE=FNAME_TDF)

DO i = 1, NMIN
    IF (MIN_NODES(i)%INFASTESTCURR .AND. (.NOT. DOALL)) THEN
        WRITE(7,*) MIN_NODES(i)%MIN_ID
        MIN_NODES(i)%INFASTESTCURR = .FALSE.
    ELSE IF (MIN_NODES(i)%INFASTESTALL .AND. DOALL) THEN
        WRITE(7,*) MIN_NODES(i)%MIN_ID
    ENDIF
ENDDO
DO i = 1, 2*NTS
    IF (TS_EDGES(i)%INFASTESTCURR .AND. (.NOT. DOALL)) THEN
        WRITE(8,*) TS_EDGES(i)%TS_ID
        TS_EDGES(i)%INFASTESTCURR = .FALSE.
    ELSE IF (TS_EDGES(i)%INFASTESTALL .AND. DOALL) THEN
        WRITE(8,*) TS_EDGES(i)%TS_ID
    ENDIF
ENDDO

CLOSE(7)
CLOSE(8)

END SUBROUTINE PRINT_INFASTEST

END MODULE KDP_SUBS

!-------------------------------------------------------------------------------
! Subroutine for finding the k shortest "distinct" paths
!-------------------------------------------------------------------------------
SUBROUTINE KDISTINCTPATHS()

USE PORFUNCS
USE COMMONS
USE GRAPH_KDP
USE TREE_KDP
USE PRIORITY_QUEUE_KDP
USE KDP_SUBS
IMPLICIT NONE
EXTERNAL MAKED4
EXTERNAL RATECONST_SETUP

INTEGER :: i, j, k, Z, TS_ID
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

PRINT *, 'kdistinctpaths> doing KDISTINCTPATHS'
PRINT *, 'finding the ', KPATHS, ' best paths'
PRINT *, 'there are ', NMIN, ' minima and ', NTS, ' transition states'
IF (KPATHS < 2) THEN
    PRINT *, 'kdistinctpaths> error: KDISTINCTPATHS keyword is used to find two or more paths'
    STOP
ENDIF

LJ1 = LOCATIONA(1)
LJ2 = LOCATIONB(1)
NOMOREPATHS = .FALSE.
NULLIFY(FIRST_TS_PATH)
NULLIFY(A_TS_PATH)
NULLIFY(HEADLLRED,DYNLLRED)

PRINT *, 'reactant min: ', LJ1, 'product min: ', LJ2

! allocate the arrays of derived data types
IF (ALLOCATED(MIN_NODES)) DEALLOCATE(MIN_NODES)
ALLOCATE(MIN_NODES(NMIN))
IF (ALLOCATED(TS_EDGES)) DEALLOCATE(TS_EDGES)
ALLOCATE(TS_EDGES(2*NTS))
IF (ALLOCATED(SP_TREE)) DEALLOCATE(SP_TREE)
ALLOCATE(SP_TREE(NMIN))

!!! build the array DMATMC containing TS weights
PRINT *, 'kdistinctpaths> setting up rate constants'
CALL RATECONST_SETUP(KSUM,DEADTS,NDEAD,.FALSE.,CUT_UNDERFLOW)
PRINT *, 'kdistinctpaths> setting up adjacency matrix'
CALL MAKED4(DMATMC,NCOL,NVAL,DEADTS,KSUM,INDEX_TS)
PRINT *, 'kdistinctpaths> number of dead TSs: ', NDEAD


!!! translate input minima and transition state arrays into the GRAPH_KDP derived type

DO i = 1, NMIN ! Nullify all pointers for nodes (representing minima)
    NULLIFY(MIN_NODES(i)%TOP_TO,MIN_NODES(i)%TOP_FROM)!,MIN_NODES(i)%NEXT_TO,MIN_NODES(i)%NEXT_FROM)
    NULLIFY(SP_TREE(i)%PARENT_NODE,SP_TREE(i)%PARENT_TS)
ENDDO

DO i = 1, 2*NTS ! Nullify all pointers for edges (representing transition states) (NB bidirectional)
    NULLIFY(TS_EDGES(i)%FROM_NODE,TS_EDGES(i)%TO_NODE,TS_EDGES(i)%NEXT_TO,TS_EDGES(i)%NEXT_FROM)
    TS_EDGES(i)%CHILD = .FALSE.
    TS_EDGES(i)%INFASTESTALL = .FALSE.
    TS_EDGES(i)%INFASTESTCURR = .FALSE.
ENDDO

! node IDs are equal to positions of minima in the array read from min.data
PRINT *, 'kdistinctpaths> assigning node information to array'
DO i = 1, NMIN
    MIN_NODES(i)%MIN_ID = 0
    MIN_NODES(i)%VISITED = .FALSE.
    MIN_NODES(i)%RED = .FALSE.
    MIN_NODES(i)%INFASTESTALL = .FALSE.
    MIN_NODES(i)%INFASTESTCURR = .FALSE.
    IF (NCOL(i)==0) CYCLE ! min has no neighbours - MIN_ID = 0 is therefore used to signal an invalid node
    MIN_NODES(i)%MIN_ID = i
ENDDO

PRINT *, 'kdistinctpaths> assigning edge information to array'
! edge IDs are equal to positions of minima in the array read from ts.data
DO i = 1, NTS
    TS_EDGES(i)%TS_ID = 0
    TS_EDGES(i+NTS)%TS_ID = 0
    IF (DEADTS(i)) CYCLE
    TS_EDGES(i)%TS_ID = i
    TS_EDGES(i+NTS)%TS_ID = i
    ! assign edge weights and to/from nodes in derived data type
    IF ((PLUS(i)==LJ1) .OR. (MINUS(i)==LJ2)) THEN
        TS_EDGES(i)%W = -LOG(DMATMC(i,1))
        TS_EDGES(i)%FROM_NODE => MIN_NODES(PLUS(i))
        TS_EDGES(i)%TO_NODE => MIN_NODES(MINUS(i))
    ELSEIF ((MINUS(i)==LJ1) .OR. (PLUS(i)==LJ2)) THEN
        TS_EDGES(NTS+i)%W = -LOG(DMATMC(i,2))
        TS_EDGES(NTS+i)%FROM_NODE => MIN_NODES(MINUS(i))
        TS_EDGES(NTS+i)%TO_NODE => MIN_NODES(PLUS(i))
    ElSE
        TS_EDGES(i)%W = -LOG(DMATMC(i,1))
        TS_EDGES(i)%FROM_NODE => MIN_NODES(PLUS(i))
        TS_EDGES(i)%TO_NODE => MIN_NODES(MINUS(i))

        TS_EDGES(NTS+i)%W = -LOG(DMATMC(i,2))
        TS_EDGES(NTS+i)%FROM_NODE => MIN_NODES(MINUS(i))
        TS_EDGES(NTS+i)%TO_NODE => MIN_NODES(PLUS(i))
    ENDIF
ENDDO

IF (KDPDUMPEDGEST) THEN ! dump the edge weights to a file kdp_tsedges.dat to be read by Python script
    OPEN(7,FILE="kdp_tsedges.dat")
    DO i = 1, NTS
        WRITE(7,*) TS_EDGES(i)%W, TS_EDGES(i+NTS)%W ! NB zero weights indicate "dead" TS
    ENDDO
    CLOSE(7)
    PRINT *, 'kdistinctpaths> finished printing edge weights to file as requested, stopping'
    STOP
ENDIF

PRINT *, 'kdistinctpaths> building graph data structure to represent the KTN'
DO i = 1, NMIN
    IF (NCOL(i) == 0) CYCLE ! min has no neighbours
    DO j = 1, NCOL(i) ! loop over all neighbours of min i
        IF (ASSOCIATED(TS_EDGES(INDEX_TS(j,i))%TO_NODE)) THEN
            IF (TS_EDGES(INDEX_TS(j,i))%TO_NODE%MIN_ID == i) THEN
                CALL ADD_TO_EDGE(i, INDEX_TS(j,i))
            ELSEIF (TS_EDGES(INDEX_TS(j,i))%FROM_NODE%MIN_ID == i) THEN
                CALL ADD_FROM_EDGE(i, INDEX_TS(j,i))
            ELSE
                PRINT *, 'kdistinctpaths> something went wrong'
                STOP
            ENDIF
        ENDIF
        IF (ASSOCIATED(TS_EDGES(INDEX_TS(j,i)+NTS)%TO_NODE)) THEN
            IF (TS_EDGES(INDEX_TS(j,i)+NTS)%TO_NODE%MIN_ID == i) THEN
               CALL ADD_TO_EDGE(i, INDEX_TS(j,i)+NTS)
            ELSEIF (TS_EDGES(INDEX_TS(j,i) + NTS)%FROM_NODE%MIN_ID == i) THEN
               CALL ADD_FROM_EDGE(i, INDEX_TS(j,i)+NTS)
            ELSE
               PRINT *, 'kdistinctpaths> something went wrong'
               STOP
            ENDIF
        ENDIF
    ENDDO
ENDDO
PRINT *, 'kdistinctpaths> finished building graph data structure'

IF (DEBUG) THEN
    PRINT *, 'kdistinctpaths> running debug tests to verify implementation of data structures'
    ! simple tests node data structure is correct
    PRINT *, 'ID of node 5000: ', MIN_NODES(5000)%MIN_ID, ' energy of node 5000: ', EMIN(MIN_NODES(5000)%MIN_ID)
    ! simple test edge data structure is correct
    PRINT *, 'Is TS #5000 dead? ', DEADTS(5000)
    PRINT *, 'ID of TS 6000: ', TS_EDGES(6000)%TS_ID, ' energy of TS 6000: ', ETS(TS_EDGES(6000)%TS_ID), ' weight: ', &
             TS_EDGES(6000)%W, 'from node: ', TS_EDGES(6000)%FROM_NODE%MIN_ID, &
             ' to node: ', TS_EDGES(6000)%TO_NODE%MIN_ID
    
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
    PRINT *, 'Now I am going to delete the top FROM edge from node',i,'in the GRAPH data structure: '
    CALL DEL_FROM_EDGE(i)
    PRINT *, 'The FROM edges for node',i,'are now:'
    TSEDGEPTR => MIN_NODES(i)%TOP_FROM
    DO WHILE (ASSOCIATED(TSEDGEPTR))
        PRINT *, 'TS ID: ', TSEDGEPTR%TS_ID, 'from ', TSEDGEPTR%FROM_NODE%MIN_ID, 'to ', TSEDGEPTR%TO_NODE%MIN_ID
        TSEDGEPTR => TSEDGEPTR%NEXT_FROM
    ENDDO
    NULLIFY(TSEDGEPTR)
    
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
ENDIF

!!! Now the algorithm of Frigioni, Marchetti-Spaccamella & Nanni
! First find the initial shortest path tree by Dijkstra's algorithm
CALL DIJKSTRA_KDP(LJ1)

MIN_NODES(LJ2)%INFASTESTALL = .TRUE.
NULLIFY(RLEDGEPTR)
CALL PRINT_PATH_KDP(RLEDGEPTR,LJ2,1)
CALL PRINT_INFASTEST(LJ2,1)

DO k = 2, KPATHS
    ! block the rate-limiting edge in both directions
    TS_EDGES(RLEDGEPTR%TS_ID)%W = HUGE(0.0D0)
    TS_EDGES(RLEDGEPTR%TS_ID+NTS)%W = HUGE(0.0D0)
    ! queue the owner node of the rate-limiting edge
    CALL KDP_PQ%ENQUEUE(SP_TREE(RLEDGEPTR%TO_NODE%MIN_ID)%CURR_DIST,RLEDGEPTR%TO_NODE)
    NULLIFY(RLEDGEPTR)

    ! Now find nodes in the shortest path tree where we need to check for alternative connections
    CALL MARCHETTI_COLOURING()
    CALL KDP_PQ%DESTROY()
    CALL MARCHETTI_PROCESS_REDS(LJ2,LJ1)
    CALL KDP_PQ%DESTROY()

    CALL PRINT_PATH_KDP(RLEDGEPTR,LJ2,k)
    CALL PRINT_INFASTEST(LJ2,k)

    DO i = 1, NMIN
        MIN_NODES(i)%RED = .FALSE.
    ENDDO
ENDDO
NULLIFY(RLEDGEPTR)

CALL PRINT_INFASTEST(LJ2,0,.TRUE.)

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
IF (ALLOCATED(MIN_NODES)) DEALLOCATE(MIN_NODES)
IF (ALLOCATED(TS_EDGES)) DEALLOCATE(TS_EDGES)
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

INTEGER :: i, n
LOGICAL :: NONRED_SHORTER
TYPE(PQ_ENTRY) :: PQENTRYOBJ
TYPE(NODE), POINTER :: MINNODEPTR
TYPE(EDGE), POINTER :: TSEDGEPTR

NULLIFY(MINNODEPTR,TSEDGEPTR)
NULLIFY(HEADLLRED,DYNLLRED)

n = 0
DO WHILE (KDP_PQ%PQ_SZ>0)
    NONRED_SHORTER = .FALSE.
    PQENTRYOBJ = KDP_PQ%TOP()
    MINNODEPTR => PQENTRYOBJ%PQ_NODEPTR
    TSEDGEPTR => MINNODEPTR%TOP_FROM
    ! looping over FROM edges to find any NONRED_SHORTER nodes
    DO WHILE (ASSOCIATED(TSEDGEPTR))
        IF (.NOT. MINNODEPTR%RED .AND. TSEDGEPTR%W>1.D-10 .AND. & ! quack check this next statement
            SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%CURR_DIST+TSEDGEPTR%W==SP_TREE(TSEDGEPTR%FROM_NODE%MIN_ID)%CURR_DIST) THEN
            NONRED_SHORTER = .TRUE.
            ! here vertex should be coloured 'pink', but with double precision edges this situation rarely happens
        ENDIF
        MINNODEPTR => TSEDGEPTR%FROM_NODE
        TSEDGEPTR => TSEDGEPTR%NEXT_FROM
    ENDDO
    NULLIFY(MINNODEPTR,TSEDGEPTR)
    IF (.NOT. NONRED_SHORTER) THEN
        MINNODEPTR => PQENTRYOBJ%PQ_NODEPTR
        IF (.NOT. MINNODEPTR%RED) THEN
            n = n + 1
        ELSE
            PRINT *, 'kdistinctpaths> something went wrong: min',MINNODEPTR%MIN_ID,'is already red'
            STOP
        ENDIF
        MINNODEPTR%RED = .TRUE.
        ! keep a linked list of the nodes coloured red
        IF (.NOT. ASSOCIATED(HEADLLRED)) THEN
            ALLOCATE(HEADLLRED)
            NULLIFY(HEADLLRED%NODE_INLIST,HEADLLRED%NEXTLLENTRY)
            HEADLLRED%NODE_INLIST => MINNODEPTR
            ALLOCATE(HEADLLRED%NEXTLLENTRY)
            DYNLLRED => HEADLLRED%NEXTLLENTRY
        ELSE
            DYNLLRED%NODE_INLIST => MINNODEPTR
            ALLOCATE(DYNLLRED%NEXTLLENTRY)
            DYNLLRED => DYNLLRED%NEXTLLENTRY
        ENDIF
        TSEDGEPTR => MINNODEPTR%TOP_FROM
        ! trace the shortest path tree to find the children of z and add to the min priority queue
        DO WHILE (ASSOCIATED(TSEDGEPTR))
            IF (TSEDGEPTR%CHILD .AND. .NOT. TSEDGEPTR%W>=HUGE(0.0D0)) THEN ! edge is in the shortest path tree
                MINNODEPTR => TSEDGEPTR%TO_NODE
                CALL KDP_PQ%ENQUEUE(SP_TREE(MINNODEPTR%MIN_ID)%CURR_DIST,MINNODEPTR)
            ENDIF
            TSEDGEPTR => TSEDGEPTR%NEXT_FROM
        ENDDO
    ENDIF
    NULLIFY(DYNLLRED%NEXTLLENTRY)
ENDDO

!PRINT *, 'length of linked list (number of RED nodes):', n
!PRINT *, 'the child nodes are:'
!NULLIFY(DYNLLRED)
!DYNLLRED => HEADLLRED
!DO WHILE (ASSOCIATED(DYNLLRED%NEXTLLENTRY))
!    PRINT *, DYNLLRED%NODE_INLIST%MIN_ID
!    DYNLLRED => DYNLLRED%NEXTLLENTRY
!ENDDO

NULLIFY(DYNLLRED)
NULLIFY(MINNODEPTR,TSEDGEPTR)

END SUBROUTINE MARCHETTI_COLOURING


!-------------------------------------------------------------------------------
! Subroutine for second step of Marchetti-Spaccamella algorithm:
! Processes vertices that are coloured "red" to find best alternative parents and
! make updates to the shortest path tree accordingly
!-------------------------------------------------------------------------------
SUBROUTINE MARCHETTI_PROCESS_REDS(LJ2,LJ1)

USE COMMONS
USE GRAPH_KDP
USE TREE_KDP
USE PRIORITY_QUEUE_KDP
IMPLICIT NONE

INTEGER :: LJ2, LJ1, n
LOGICAL :: NONRED_NBR
DOUBLE PRECISION :: FLVLZ, FLVLZBEST, EDGECOST
TYPE(PQ_ENTRY) :: PQENTRYOBJ
TYPE(NODE), POINTER :: MINNODEPTR, BESTNONREDNBR
TYPE(EDGE), POINTER :: TSEDGEPTR, BESTNBREDGE

n = 0
NULLIFY(MINNODEPTR,TSEDGEPTR)
DYNLLRED => HEADLLRED
DO WHILE (ASSOCIATED(DYNLLRED%NEXTLLENTRY))
    NULLIFY(BESTNONREDNBR,BESTNBREDGE)
    FLVLZBEST = HUGE(0.D0)
    NONRED_NBR = .FALSE.
    MINNODEPTR => DYNLLRED%NODE_INLIST
    ! loop over neighbours of the current red node
    TSEDGEPTR => MINNODEPTR%TOP_TO
    DO WHILE (ASSOCIATED(TSEDGEPTR))
        IF (.NOT. TSEDGEPTR%FROM_NODE%RED) THEN
            NONRED_NBR = .TRUE.
            FLVLZ = SP_TREE(TSEDGEPTR%FROM_NODE%MIN_ID)%CURR_DIST + TSEDGEPTR%W
            IF (.NOT. ASSOCIATED(BESTNONREDNBR) .OR. FLVLZ < FLVLZBEST) THEN
                BESTNONREDNBR => TSEDGEPTR%FROM_NODE
                BESTNBREDGE => TSEDGEPTR
                EDGECOST = TSEDGEPTR%W
                FLVLZBEST = FLVLZ
            ENDIF
        ENDIF
        TSEDGEPTR => TSEDGEPTR%NEXT_TO
    ENDDO
    IF (.NOT. NONRED_NBR) THEN ! disclude this node from the shortest path tree
        SP_TREE(MINNODEPTR%MIN_ID)%CURR_DIST = HUGE(0.D0)
        NULLIFY(SP_TREE(MINNODEPTR%MIN_ID)%PARENT_NODE)
    ELSE
        SP_TREE(MINNODEPTR%MIN_ID)%CURR_DIST = SP_TREE(BESTNONREDNBR%MIN_ID)%CURR_DIST + EDGECOST
        SP_TREE(MINNODEPTR%MIN_ID)%PARENT_NODE => BESTNONREDNBR
        SP_TREE(MINNODEPTR%MIN_ID)%PARENT_TS%CHILD = .FALSE. ! need to give old CHILD edge .FALSE. label
        BESTNBREDGE%CHILD = .TRUE.
        SP_TREE(MINNODEPTR%MIN_ID)%PARENT_TS => BESTNBREDGE
        CALL KDP_PQ%ENQUEUE(SP_TREE(MINNODEPTR%MIN_ID)%CURR_DIST,MINNODEPTR)
        MINNODEPTR%VISITED = .TRUE. ! flag indicates node has been queued
    ENDIF    

    ! go to next node in linked list
    n = n+1
    DYNLLRED => DYNLLRED%NEXTLLENTRY
    DEALLOCATE(HEADLLRED)
    HEADLLRED => DYNLLRED
ENDDO

NULLIFY(MINNODEPTR,TSEDGEPTR,BESTNONREDNBR,BESTNBREDGE)
DO WHILE (KDP_PQ%PQ_SZ>0)
    PQENTRYOBJ = KDP_PQ%TOP()
    MINNODEPTR => PQENTRYOBJ%PQ_NODEPTR
    IF (.NOT. MINNODEPTR%VISITED) CYCLE ! found "out-of-date" entry in priority queue
    MINNODEPTR%VISITED = .FALSE. ! flag indicates node has been dequeued
    TSEDGEPTR => MINNODEPTR%TOP_FROM
    DO WHILE (ASSOCIATED(TSEDGEPTR))
        IF (TSEDGEPTR%TO_NODE%RED) THEN
            FLVLZ = SP_TREE(MINNODEPTR%MIN_ID)%CURR_DIST + TSEDGEPTR%W
            IF (FLVLZ < SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%CURR_DIST) THEN
                SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%CURR_DIST = FLVLZ
                SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%PARENT_NODE => MINNODEPTR
                SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%PARENT_TS%CHILD = .FALSE. ! old edge in sp tree
                TSEDGEPTR%CHILD = .TRUE.
                SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%PARENT_TS => TSEDGEPTR
                ! add improved entry to priority queue
                CALL KDP_PQ%ENQUEUE(SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%CURR_DIST,TSEDGEPTR%TO_NODE)
                TSEDGEPTR%TO_NODE%VISITED = .TRUE.
            ENDIF
        ENDIF
        TSEDGEPTR => TSEDGEPTR%NEXT_FROM
    ENDDO
ENDDO

NULLIFY(MINNODEPTR)
IF (.NOT. ASSOCIATED(SP_TREE(LJ2)%PARENT_NODE)) THEN
    PRINT *, 'kdistinctpaths> no existing path to endpoint node',LJ2
    STOP
ENDIF
MINNODEPTR => SP_TREE(LJ2)%PARENT_NODE
DO WHILE (.TRUE.)
    IF (.NOT. ASSOCIATED(SP_TREE(MINNODEPTR%MIN_ID)%PARENT_NODE)) THEN
        IF (MINNODEPTR%MIN_ID /= LJ1) THEN
            PRINT *, 'kdistinctpaths> no existing path to endpoint node',LJ2
            STOP
        ELSE
            EXIT
        ENDIF
    ENDIF
    MINNODEPTR => SP_TREE(MINNODEPTR%MIN_ID)%PARENT_NODE
ENDDO

NULLIFY(MINNODEPTR,TSEDGEPTR)

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

! initialisation
! node IDs for tree data type are equal to positions of minima in the array read from min.data
DO i = 1, NMIN
    SP_TREE(i)%MIN_ID = i
    NULLIFY(SP_TREE(i)%PARENT_NODE,SP_TREE(i)%PARENT_TS)
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
    TSEDGEPTR => MIN_NODES(n)%TOP_FROM
    MIN_NODES(n)%VISITED = .TRUE.
    DO WHILE (ASSOCIATED(TSEDGEPTR)) ! loop over all neighbours of min n (FROM n TO x)
        !PRINT *, 'current edge is FROM/TO: ', TSEDGEPTR%FROM_NODE%MIN_ID,TSEDGEPTR%TO_NODE%MIN_ID
        ALT = SP_TREE(MIN_NODES(n)%MIN_ID)%CURR_DIST + TSEDGEPTR%W
        IF (ALT < SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%CURR_DIST) THEN
            SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%CURR_DIST = ALT
            SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%PARENT_NODE => MIN_NODES(n)
            IF (ASSOCIATED(SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%PARENT_TS)) THEN
                SP_TREE(TSEDGEPTR%TO_NODE%MIN_ID)%PARENT_TS%CHILD = .FALSE.
            ENDIF
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
        ENDIF
    ENDDO
ENDDO

DO i = 1, NMIN
    MIN_NODES(i)%VISITED = .FALSE.
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
