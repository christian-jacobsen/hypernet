import numpy as np

from hypernet.src.general import const
from hypernet.src.general import utils


class SpecieCoeffs(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        *args,
        **kwargs
    ):


    O2(1)+O=3O
    O2(2)+O=O2(1)+O

        for index, row in self.arrhenius_coeff.iterrows():
            rate = self.modified_arrhenius(T, row['A'], row['n'], row['Ea'])
            reac_name, *bins_idx = index.split('_')
            if reac_name == 'diss':
                idx = int(bins_idx[0])-1
            else:
                idx = int(bins_idx[0])-1, int(bins_idx[1])-1
            rates_coeff[reac_name][idx] = rate


    def dkfdT(self, T):
        return self.reactionRate.dkdT(T)


    def setLRhs(self, reacStr):
        lhs, rhs = reacStr.split('=')
        lhs, rhs = setXhs(lhs)
        lhs, rhs = setXhs(rhs)

        return self.reactionRate.dkdT(T)

    def setXhs(self, xhsStr):
        identifier = xhsStr.split('+')
        chars = iter([c for c in identifier])
        if next(chars).isnumeric():


        if (t.isNumber())
        {
            stoichCoeff = t.number();
        }
        else
        {
            stoichCoeff = 1;
        }

        return self.reactionRate.dkdT(T)

    def setRhs(self, T):
        return self.reactionRate.dkdT(T)


import re
s = '111A222B333C'
res = re.split('(\d+)', s)
print(res)
# ['', '111', 'A', '222', 'B', '333', ' C']


    Foam::specieCoeffs::specieCoeffs
    (
        const speciesTable& species,
        Istream& is
    )
    {
        token t(is);
        if (t.isNumber())
        {
            stoichCoeff = t.number();
            is >> t;
        }
        else
        {
            stoichCoeff = 1;
        }

        exponent = stoichCoeff;

        if (t.isWord())
        {
            word specieName = t.wordToken();

            size_t i = specieName.find('^');

            if (i != word::npos)
            {
                string exponentStr = specieName
                (
                    i + 1,
                    specieName.size() - i - 1
                );
                exponent = atof(exponentStr.c_str());
                specieName = specieName(0, i);
            }

            if (species.found(specieName))
            {
                index = species[specieName];
            }
            else
            {
                FatalIOErrorInFunction(is)
                    << "Specie " << specieName
                    << " not found in table " << species
                    << exit(FatalIOError);

                index = -1;
            }
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Expected a word but found " << t.info()
                << exit(FatalIOError);
        }
    }


    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    void Foam::specieCoeffs::setLRhs
    (
        Istream& is,
        const speciesTable& species,
        List<specieCoeffs>& lhs,
        List<specieCoeffs>& rhs
    )
    {
        DynamicList<specieCoeffs> dlrhs;

        while (is.good())
        {
            dlrhs.append(specieCoeffs(species, is));

            if (dlrhs.last().index != -1)
            {
                token t(is);
                if (t.isPunctuation())
                {
                    if (t == token::ADD)
                    {
                    }
                    else if (t == token::ASSIGN)
                    {
                        lhs = dlrhs.shrink();
                        dlrhs.clear();
                    }
                    else
                    {
                        rhs = dlrhs.shrink();
                        is.putBack(t);
                        return;
                    }
                }
                else
                {
                    rhs = dlrhs.shrink();
                    is.putBack(t);
                    return;
                }
            }
            else
            {
                dlrhs.remove();
                if (is.good())
                {
                    token t(is);
                    if (t.isPunctuation())
                    {
                        if (t == token::ADD)
                        {
                        }
                        else if (t == token::ASSIGN)
                        {
                            lhs = dlrhs.shrink();
                            dlrhs.clear();
                        }
                        else
                        {
                            rhs = dlrhs.shrink();
                            is.putBack(t);
                            return;
                        }
                    }
                }
                else
                {
                    if (!dlrhs.empty())
                    {
                        rhs = dlrhs.shrink();
                    }
                    return;
                }
            }
        }

        FatalIOErrorInFunction(is)
            << "Cannot continue reading reaction data from stream"
            << exit(FatalIOError);
    }


    void Foam::specieCoeffs::reactionStr
    (
        OStringStream& reaction,
        const speciesTable& species,
        const List<specieCoeffs>& scs
    )
    {
        for (label i = 0; i < scs.size(); ++i)
        {
            if (i > 0)
            {
                reaction << " + ";
            }
            if (mag(scs[i].stoichCoeff - 1) > small)
            {
                reaction << scs[i].stoichCoeff;
            }
            reaction << species[scs[i].index];
            if (mag(scs[i].exponent - scs[i].stoichCoeff) > small)
            {
                reaction << "^" << scs[i].exponent;
            }
        }
    }


    Foam::string Foam::specieCoeffs::reactionStr
    (
        OStringStream& reaction,
        const speciesTable& species,
        const List<specieCoeffs>& lhs,
        const List<specieCoeffs>& rhs
    )
    {
        reactionStr(reaction, species, lhs);
        reaction << " = ";
        reactionStr(reaction, species, rhs);
        return reaction.str();
    }


    // ************************************************************************* //


