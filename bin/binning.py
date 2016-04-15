

def CountBin(run_ID, ML_bin, SA_bin, VF_bin):

    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker
    from runDB_declarative import RunData, Base

    engine = create_engine( "sqlite:///HTSOHM-dev.db" )
    Base.metadata.bind = engine

    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    BinCount = session.query(RunData).filter(RunData.Run == run_ID,
                                             RunData.Bin_ML == ML_bin,
                                             RunData.Bin_SA == SA_bin,
                                             RunData.Bin_VF == VF_bin).count()

    return BinCount


def CountAll(run_ID):
    
    import numpy as np

    bins = int(run_ID[-1])
    AllCounts = np.zeros([bins, bins, bins])

    for i in range(bins):
        for j in range(bins):
            for k in range(bins):
                
                b_count = CountBin(run_ID, i, j, k)
                AllCounts[i,j,k] = b_count

    return AllCounts


def SelectParents(run_ID, children_per_generation, generation):

    import numpy as np

    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker
    from runDB_declarative import RunData, Base

    engine = create_engine( "sqlite:///HTSOHM-dev.db" )
    Base.metadata.bind = engine

    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    bins = int(run_ID[-1])
    counts = CountAll(run_ID)
    weights = np.zeros([bins, bins, bins])

    for i in range(bins):
        for j in range(bins):
            for k in range(bins):    
                if counts[i,j,k] != 0.:
                    weights[i,j,k] = counts.sum() / counts[i,j,k]
    weights = weights / weights.sum()

    w_list = []
    ID_list = []
    for i in range(bins):
        for j in range(bins):
            w_list = np.concatenate( [w_list, weights[i,j,:]] )
            for k in range(bins):
                bin_IDs = []
                res = session.query(RunData).filter(RunData.Run == run_ID,
                                                         RunData.Bin_ML == i,
                                                         RunData.Bin_SA == j,
                                                         RunData.Bin_VF == k).all()
                for item in res:
                    bin_IDs.append(item.Mat)
                ID_list = ID_list + [bin_IDs]

    first = generation * children_per_generation
    last = (generation + 1) * children_per_generation
    new_mat_IDs = np.arange(first, last)                   # IDs for next generation of materials

    for i in new_mat_IDs:

        p_bin = np.random.choice(ID_list, p=w_list)
        p_ID = np.random.choice(p_bin)                     # Select parent for new material

        # Write parent IDs to database...
        check_first = session.query(RunData).filter( RunData.Run == run_ID,
                                                     RunData.Mat == str(i) ).count()
        if not check_first:
            new_mat = RunData( Run=run_ID, Mat=str(i) )
            session.add(new_mat)
            session.commit()

        NextMat = session.query(RunData).filter( RunData.Run == run_ID,
                                                 RunData.Mat == str(i) )
        NextMat.update({'Parent': str(p_ID)})
        session.commit()
    

def CheckConvergance(run_ID):

    counts = CountAll(run_ID)
       
