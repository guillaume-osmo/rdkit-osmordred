# Guide de compilation pour le patch calcRNCG_RPCG

## ✅ Patch appliqué

Le patch pour `calcRNCG_RPCG` a été appliqué avec succès dans :
- `Code/GraphMol/Descriptors/Osmordred.cpp` (ligne ~7457)

## 🔨 Compilation

### Option 1: Compilation complète (si RDKit n'est pas déjà compilé)

```bash
cd /Users/guillaume-osmo/Github/rdkit-pypi

# Activer l'environnement de build
source /Users/guillaume-osmo/miniconda3/etc/profile.d/conda.sh
conda activate rdkit_build_py311_complete

# Vérifier que RDKit source est à jour
cd rdkit-source
git status  # Vérifier que Osmordred.cpp est modifié

# Retourner au répertoire principal
cd ..

# Lancer le build complet
export CMAKE_BUILD_PARALLEL_LEVEL=12
bash build_local.sh
```

**Temps estimé:** 30-60 minutes

### Option 2: Bypass rebuild (si RDKit est déjà compilé)

Si RDKit a déjà été compilé précédemment, on peut bypasser la recompilation complète :

```bash
cd /Users/guillaume-osmo/Github/rdkit-pypi

# Activer l'environnement de build
source /Users/guillaume-osmo/miniconda3/etc/profile.d/conda.sh
conda activate rdkit_build_py311_complete

# Vérifier si RDKit est déjà compilé
LIB_COUNT=$(find build/temp.macosx-11.0-arm64-cpython-311/rdkit_install/lib -name "*.dylib" 2>/dev/null | wc -l | tr -d ' ')

if [ "$LIB_COUNT" -gt 200 ]; then
    echo "✅ RDKit déjà compilé ($LIB_COUNT bibliothèques) - recompilation partielle seulement"
    
    # Forcer la recompilation de Osmordred seulement
    # En supprimant les objets compilés de Descriptors
    find build/temp.macosx-11.0-arm64-cpython-311/rdkit/build/Code/GraphMol/Descriptors -name "Osmordred.o" -delete 2>/dev/null || true
    
    # Relancer le build (CMake détectera les changements et recompilera seulement ce qui est nécessaire)
    export CMAKE_BUILD_PARALLEL_LEVEL=12
    python setup.py bdist_wheel
else
    echo "❌ RDKit pas encore compilé - build complet nécessaire"
    export CMAKE_BUILD_PARALLEL_LEVEL=12
    bash build_local.sh
fi
```

**Temps estimé:** 5-10 minutes (si RDKit déjà compilé)

## 📦 Installation

Après compilation réussie :

```bash
# Activer l'environnement osmo
conda activate osmo

# Installer la nouvelle wheel
pip install --force-reinstall dist/rdkit-*cp311*.whl

# Vérifier l'installation
python -c "
from rdkit import Chem
from rdkit.Chem import Osmordred
print('✅ Osmordred importé avec succès')

# Tester avec la molécule problématique (Hg)
mol = Chem.MolFromSmiles('Cc1ccc([N+](=O)[O-])[c]2c1[O][Hg]2')
if mol:
    try:
        result = Osmordred.CalcRNCGRPCG(mol)
        print(f'✅ calcRNCG_RPCG retourné: {result}')
        import numpy as np
        if np.isnan(result[0]) and np.isnan(result[1]):
            print('✅ NaN retournés correctement pour la molécule Hg!')
        else:
            print(f'⚠️  Résultats non-NaN: {result}')
    except Exception as e:
        print(f'❌ Erreur: {e}')
else:
    print('❌ Molécule non parsée')
"
```

## 🧪 Test complet

Pour tester avec plusieurs molécules :

```bash
python3 << 'PYEOF'
from rdkit import Chem
from rdkit.Chem import Osmordred
import numpy as np

# Molécules de test
test_molecules = [
    ('CCO', 'Éthanol (normal)'),
    ('Cc1ccc([N+](=O)[O-])[c]2c1[O][Hg]2', 'Mercure (problématique)'),
    ('CC', 'Éthane (normal)'),
]

print("="*80)
print("TEST: calcRNCG_RPCG avec patch Gasteiger NaN")
print("="*80)

for smiles, desc in test_molecules:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        try:
            result = Osmordred.CalcRNCGRPCG(mol)
            has_nan = np.isnan(result[0]) or np.isnan(result[1])
            status = "✅ NaN" if has_nan else "✅ OK"
            print(f"\n{desc} ({smiles}):")
            print(f"  {status} Résultat: {result}")
        except Exception as e:
            print(f"\n{desc} ({smiles}):")
            print(f"  ❌ Erreur: {e}")

print("\n" + "="*80)
PYEOF
```

## 📝 Notes

1. **Recompilation partielle**: CMake détecte automatiquement les changements dans `Osmordred.cpp` et recompile seulement les fichiers nécessaires.

2. **Vérification du patch**: Le patch ajoute :
   - Try-catch autour de `computeGasteigerCharges`
   - Vérification des NaN dans les charges
   - Retour de `{NaN, NaN}` si Gasteiger échoue

3. **Includes**: Les includes nécessaires (`<limits>` et `<cmath>`) sont déjà présents dans le fichier.

## 🔍 Vérification du patch

Pour vérifier que le patch est bien appliqué :

```bash
cd /Users/guillaume-osmo/Github/rdkit-pypi/rdkit-source
grep -A 5 "gasteiger_failed = false" Code/GraphMol/Descriptors/Osmordred.cpp
```

Vous devriez voir le code du patch avec `try-catch` et la gestion des NaN.

