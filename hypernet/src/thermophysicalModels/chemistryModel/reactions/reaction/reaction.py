import numpy as np

from hypernet.src.general import const
from hypernet.src.general import utils


class Reaction(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        *args,
        **kwargs
    ):



    #include "reaction.H"
    #include "dictionary.H"


    // * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

    Foam::label Foam::reaction::nUnNamedReactions(0);


    // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

    Foam::label Foam::reaction::getNewReactionID()
    {
        return nUnNamedReactions++;
    }


    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    Foam::reaction::reaction
    (
        const speciesTable& species,
        const List<specieCoeffs>& lhs,
        const List<specieCoeffs>& rhs
    )
    :
        name_("un-named-reaction-" + Foam::name(getNewReactionID())),
        species_(species),
        lhs_(lhs),
        rhs_(rhs)
    {}


    Foam::reaction::reaction
    (
        const reaction& r,
        const speciesTable& species
    )
    :
        name_(r.name() + "Copy"),
        species_(species),
        lhs_(r.lhs_),
        rhs_(r.rhs_)
    {}


    Foam::reaction::reaction
    (
        const speciesTable& species,
        const dictionary& dict
    )
    :
        name_(dict.dictName()),
        species_(species)
    {
        specieCoeffs::setLRhs
        (
            IStringStream(dict.lookup("reaction"))(),
            species_,
            lhs_,
            rhs_
        );
    }


    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void Foam::reaction::write(Ostream& os) const
    {
        OStringStream reaction;
        writeEntry
        (
            os,
            "reaction",
            specieCoeffs::reactionStr(reaction, species_, lhs_, rhs_)
        );
    }


    // ************************************************************************* //
