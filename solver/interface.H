#ifndef Interface_H
#define Interface_H
namespace Foam
{
    struct IMixture
    {
        virtual scalar GetA(const scalar MolWeight, const scalar T) const = 0;
        virtual scalar Cv(const scalar P, const scalar T, const scalar MolWeight) const = 0;
        virtual scalar GetTmax(const scalar MolWeight, const scalar Uthermodinamical) const = 0;
    };
};
#endif //Interface_H