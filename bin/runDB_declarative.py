import os
import sys
from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine

Base = declarative_base()

class RunData(Base):
    __tablename__ = 'RunData'
    # Defining columns for run's data output.
    id = Column(Integer, primary_key=True)
    Run = Column(String(50))
    Mat = Column(Integer)
    Abs_cccc = Column(Float)
    Abs_ccgr = Column(Float)
    Abs_mokg = Column(Float)
    Exc_cccc = Column(Float)
    Exc_ccgr = Column(Float)
    Exc_mokg = Column(Float)
    Des_AdAd = Column(Float)
    Des_HoAd = Column(Float)
    Des_tota = Column(Float)
    SA_A2 = Column(Float)
    SA_m2cc = Column(Float)
    SA_m2gr = Column(Float)
    VF_wido = Column(Float)
    Parent = Column(Integer)
    Bin_ML = Column(Integer)
    Bin_SA = Column(Integer)
    Bin_VF = Column(Integer)
    D_pass = Column(String(1))

# Create engine to storedata in local directory .db-file.
engine = create_engine("sqlite:///HTSOHM-dev.db")

# Create table in the engine.
Base.metadata.create_all(engine)

