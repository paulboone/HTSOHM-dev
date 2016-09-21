
class Base(object):

    def clone(self):
        copy = self.__class__()
        for col in self.__table__.columns:
            val = getattr(self, col.name)
            setattr(copy, col.name, val)
        return copy

from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base(cls=Base)
