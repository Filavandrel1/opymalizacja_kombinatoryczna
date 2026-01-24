#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib> 
#include <limits>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <vector>
#include <random>
#include <queue>
#include <functional>
#include <cmath>
using namespace std;


struct Parametry {
    int n = 0; 
    int k = 0;
    int delta_k = 0;
    int l_neg = 0;
    int l_poz = 0;
    bool repAllowed = true;
    bool probablePositive = false;
} P;

struct Instancja {
    string sekwencja;
    vector<string> spektrum;
    vector<string> spektrumZBledami;
    string pierwszyOligo;
    bool sekwencjaIstnieje = false;
    bool spektrumIstnieje = false;
    bool spektrumZBledamiIstnieje = false;
} I;

const char NUKLEOTYDY[4] = { 'A', 'C', 'G', 'T' };

void generujSpektrumBezBledow();
void zastosujBledyWSpektrum();
int policzNajdluzszyOverlapOgraniczony(const string& aktualnaSekwencja, const string& kandydat);
string zlozSekwencjeNaiwnieZachlannie(const Instancja& instancja);
int policzOdlegloscLevenshteina(const string& sekwencjaPrawdziwa, const string& sekwencjaOdtworzona);
pair<int, double> policzPokrycieSpektrum(const vector<string>& spektrumZBledami, const string& sekwencjaOdtworzona);
void ustawParametryDomyslneDlaAktualnejSekwencji();

// (brak dodatkowych forward declarations; tryb batch jest zdefiniowany na dole pliku)

int znajdzIndeksWierzcholka(const vector<string>& wierzcholki, const string& oligo) {
    for (int i = 0; i < wierzcholki.size(); i++) {
        if (wierzcholki[i] == oligo) {
            return i;
        }
    }
    return -1;
}

struct KrawedzGrafu {
    int doWierzcholka = -1;
    int overlap = 0;
    int kategoriaWagi = 3;
};

struct ParametryMrowkowe {
    int liczbaIteracji = 300;
    int liczbaMrowek = 150;
    int limitKandydatow = 40;
    double alfa = 1.5;
    double beta = 1.8;
    double parowanie = 0.10;
    double q = 3.0;

    int procentMrowekDijkstra = 10;  
    int coIleKrokowDijkstra = 50;
};

int policzNajdluzszyOverlapDowolny(const string& lewy, const string& prawy) {
    int dlugoscLewy = lewy.size();
    int dlugoscPrawy = prawy.size();
    if (dlugoscLewy == 0 || dlugoscPrawy == 0) {
        return 0;
    }

    int rzeczywistyMax = min({ dlugoscLewy, dlugoscPrawy });
    for (int overlap = rzeczywistyMax; overlap >= 1; --overlap) {
        bool pasuje = true;
        for (int i = 0; i < overlap; ++i) {
            if (lewy[dlugoscLewy - overlap + i] != prawy[i]) {
                pasuje = false;
                break;
            }
        }
        if (pasuje) {
            return overlap;
        }
    }
    return 0;
}

int policzKategorieKrawedzi(int overlap, int k, int delta_k) {
    if (overlap <= 0) {
        return 3;
    }
    int minimalnyDobry = max(1, (k - 1) - delta_k);
    if (overlap >= minimalnyDobry) {
        return 1;
    }
    return 2;
}

double policzHeurystyke(int overlap, int kategoria, bool czyOdwiedzony) {
    // Preferujemy overlap w granicach i nieodwiedzone wierzcholki.
    double bonusOdwiedzin = czyOdwiedzony ? 0.65 : 1.0;
    double karaKategorii = (kategoria == 1) ? 1.0 : (kategoria == 2 ? 2.5 : 6.0);
    return bonusOdwiedzin * (overlap + 1.0) / karaKategorii;
}

int policzOdlegloscHammingaDoLimitu(const string& a, const string& b, int limit) {
    if (a.size() != b.size()) return limit + 1;
    int roznice = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) {
            ++roznice;
            if (roznice > limit) return roznice;
        }
    }
    return roznice;
}

struct WynikSciezki {
    string sekwencja;
    vector<int> sciezkaWierzcholkow;
    vector<pair<int, int>> uzyteKrawedzie; // (z, do)
    int liczbaKrawedziKat1 = 0;
    int liczbaKrawedziKat2 = 0;
    int liczbaKrawedziKat3 = 0;
    int unikalneWierzcholki = 0;
    int powtorzeniaWierzcholkow = 0;
};

double policzOceneSciezki(const WynikSciezki& wynik) {
    // Kryteria z opisu:
    // - wiecej lukow wagi 1 lepsze
    // - wieksza czesc grafu uzyta lepsze
    // - wiecej lukow wagi 2/3 gorsze
    // - gdy brak 2/3, karzemy powtorki wierzcholkow
    double ocena = 0.0;
    // Uwaga: zbyt duza waga kat1 powoduje zapadanie sie w cykle (malo unikalnych wierzcholkow).
    // Dlatego wzmacniamy skladnik "uzyta czesc grafu" i zmniejszamy dominacje kat1.
    ocena += 150.0 * wynik.liczbaKrawedziKat1;
    ocena += 30.0 * wynik.unikalneWierzcholki;
    ocena -= 220.0 * wynik.liczbaKrawedziKat2;
    ocena -= 450.0 * wynik.liczbaKrawedziKat3;
    ocena -= 35.0 * wynik.powtorzeniaWierzcholkow;
    if (wynik.liczbaKrawedziKat2 + wynik.liczbaKrawedziKat3 == 0) {
        ocena -= 25.0 * wynik.powtorzeniaWierzcholkow;
    }
    return ocena;
}

vector<int> dijkstraSciezka(
    const vector<vector<KrawedzGrafu>>& graf,
    int start,
    int cel,
    int k,
    int delta_k
) {
    int n = graf.size();
    const int INF = 1e9;
    vector<int> dist(n, INF);
    vector<int> prev(n, -1);
    using Wezel = pair<int, int>; // (dist, v)
    priority_queue<Wezel, vector<Wezel>, greater<Wezel>> kolejka;

    dist[start] = 0;
    kolejka.push({ 0, start });

    while (!kolejka.empty()) {
        Wezel top = kolejka.top();
        int d = top.first;
        int v = top.second;
        kolejka.pop();
        if (d != dist[v]) continue;
        if (v == cel) break;

        for (const auto& e : graf[v]) {
            // Koszt: kategoria (1/2/3) oraz preferencja na dluzszy overlap (mniejszy koszt)
            int kosztKategorii = (e.kategoriaWagi == 1) ? 1 : (e.kategoriaWagi == 2 ? 3 : 6);
            int kosztOverlap = std::max(0, (k + delta_k) - e.overlap);
            int koszt = kosztKategorii * 10 + kosztOverlap;

            if (dist[e.doWierzcholka] > dist[v] + koszt) {
                dist[e.doWierzcholka] = dist[v] + koszt;
                prev[e.doWierzcholka] = v;
                kolejka.push({ dist[e.doWierzcholka], e.doWierzcholka });
            }
        }
    }

    if (prev[cel] == -1 && start != cel) {
        return {};
    }

    vector<int> sciezka;
    int v = cel;
    sciezka.push_back(v);
    while (v != start) {
        v = prev[v];
        if (v == -1) {
            return {};
        }
        sciezka.push_back(v);
    }
    reverse(sciezka.begin(), sciezka.end());
    return sciezka;
}

void zapiszWynikDoPliku(const string& nazwaPliku, const string& sekwencjaOdtworzona, int odlegloscLevenshteina, const pair<int, double>& wynikPokrycia) {
    ofstream plikWyjsciowy(nazwaPliku);
    if (!plikWyjsciowy) {
        cout << "Blad otwarcia pliku do zapisu." << endl;
        return;
    }

    plikWyjsciowy << "Oryginalne DNA:\n" << I.sekwencja << "\n\n";
    plikWyjsciowy << "Odtworzone DNA:\n" << sekwencjaOdtworzona << "\n\n";
    plikWyjsciowy << "Odleglosc Levenshteina: " << odlegloscLevenshteina << "\n";
    plikWyjsciowy << "Pokrycie spektrum: " << wynikPokrycia.first << "/" << I.spektrumZBledami.size()
                  << " (" << wynikPokrycia.second << "% )\n\n";
    plikWyjsciowy << "Parametry instancji:\n";
    plikWyjsciowy << "n = " << P.n << "\n";
    plikWyjsciowy << "k = " << P.k << "\n";
    plikWyjsciowy << "delta_k = " << P.delta_k << "\n";
    plikWyjsciowy << "l_neg = " << P.l_neg << "\n";
    plikWyjsciowy << "l_poz = " << P.l_poz << "\n";
    plikWyjsciowy << "repAllowed = " << (P.repAllowed ? 1 : 0) << "\n";
    plikWyjsciowy << "probablePositive = " << (P.probablePositive ? 1 : 0) << "\n";
}

void ustawParametryDomyslneDlaAktualnejSekwencji() {
    P.n = I.sekwencja.size();
    P.k = 8;
    P.delta_k = 2;
    P.l_neg = 0;
    P.l_poz = 0;
    P.repAllowed = true;
    P.probablePositive = false;
}

void wyswietlanieMenu() {
    if (I.sekwencjaIstnieje) {
        cout << "USTAWIENIA INSTANCJI" << endl;
        cout << "n = " << P.n
            << ", k = " << P.k
            << ", delta_k = " << P.delta_k
            << ", l_neg = " << P.l_neg
            << ", l_poz = " << P.l_poz
            << ", repAllowed = " << (P.repAllowed ? "true" : "false")
            << ", probablePositive = " << (P.probablePositive ? 1 : 0)
            << endl;
    }
    else {
        cout << "Brak wybranej sekwencji." << endl << endl;
    }
    cout << "----------MENU----------" << endl;
    cout << "1. Generator instancji" << endl;
    cout << "2. Algorytm naiwny" << endl;
    cout << "3. Metaheurystyka" << endl;
    cout << "0. Wyjscie z programu" << endl;
    cout << "------------------------" << endl;
}
void wczytajZPliku() {
    cout << "Podaj nazwe pliku do wczytania: ";
    string nazwa;
    cin >> nazwa;
    ifstream in(nazwa);
    if (!in) {
        cout << "Blad otwarcia pliku." << endl;
        return;
    }
    I.sekwencja.clear();

    string dnaLinia;
    in >> ws;
    if (!getline(in, dnaLinia)) {
        cout << "Plik pusty lub blad formatu." << endl;
        return;
    }

    if (dnaLinia.empty()) {
        cout << "Plik pusty lub blad formatu." << endl;
        return;
    }

    I.sekwencja = dnaLinia;
    I.sekwencjaIstnieje = true;

    string paramLinia;
    if (getline(in, paramLinia) && !paramLinia.empty()) {
        istringstream iss(paramLinia);
        int n, k, d, lneg, lpoz, rep, prob;
        if (iss >> n >> k >> d >> lneg >> lpoz >> rep >> prob) {
            P.n = n;
            P.k = k;
            P.delta_k = d;
            P.l_neg = lneg;
            P.l_poz = lpoz;
            P.repAllowed = (rep != 0);
            P.probablePositive = (prob != 0);
        }
        else {
            ustawParametryDomyslneDlaAktualnejSekwencji();
        }
    }
    else {
        ustawParametryDomyslneDlaAktualnejSekwencji();
    }

    // Pierwszy oligo zawsze jest prefiksem DNA; jego dlugosc nie zmienia sie w generatorze (i == 0).
    I.pierwszyOligo = (P.k > 0 && I.sekwencja.size() >= P.k)
        ? I.sekwencja.substr(0, P.k)
        : "";

    cout << "Wczytano sekwencje DNA o dlugosci " << I.sekwencja.size() << endl;
    
    generujSpektrumBezBledow();
    I.spektrumZBledami.clear();
    I.spektrumZBledamiIstnieje = false;

    int liczbaOligoZPliku = 0;
    if (in >> liczbaOligoZPliku && liczbaOligoZPliku > 0) {
        string usuniecieZnaku;
        getline(in, usuniecieZnaku);
        for (int i = 0; i < liczbaOligoZPliku; ++i) {
            string oligo;
            if (!getline(in, oligo)) break;
            if (!oligo.empty()) {
                I.spektrumZBledami.push_back(oligo);
            }
        }
        sort(I.spektrumZBledami.begin(), I.spektrumZBledami.end());
        if (!I.spektrumZBledami.empty()) {
            I.spektrumZBledamiIstnieje = true;
        }
    }

    if (!I.spektrumZBledamiIstnieje) {
        zastosujBledyWSpektrum();
    }
}

int wczytajIntZDomyslna(const string& prompt, int domyslna) {
    string linia;
    cout << prompt;
    getline(cin, linia);
    if (linia.empty()) {
        return domyslna;
    }
    return stoi(linia);
}

void pobierzParametryInstancji() {
    P.n = wczytajIntZDomyslna("Podaj n (dlugosc DNA, domyslnie 400): ", 400);
    P.k = wczytajIntZDomyslna("Podaj k (dlugosc oligo, domyslnie 8): ", 8);
    P.delta_k = wczytajIntZDomyslna("Podaj delta_k (0-2, domyslnie 2): ", 2);
    P.l_neg = wczytajIntZDomyslna("Podaj l_neg (liczba bledow negatywnych, domyslnie 0): ", 0);
    P.l_poz = wczytajIntZDomyslna("Podaj l_poz (liczba bledow pozytywnych, domyslnie 0): ", 0);

    if (P.n < 300) P.n = 300;
    if (P.n > 700) P.n = 700;

    if (P.k < 7) P.k = 7;
    if (P.k > 10) P.k = 10;

    if (P.delta_k < 0) P.delta_k = 0;
    if (P.delta_k > 2) P.delta_k = 2;

    if (P.l_neg != 0 && P.l_neg < 10) P.l_neg = 10;
    if (P.l_poz != 0 && P.l_poz < 10) P.l_poz = 10;

    string linia;
    cout << "repAllowed (dopuszczalne powtorzenia) [T/N, ENTER = domyslnie T]: ";
    getline(cin, linia);
    if (linia.empty() || linia == "T" || linia == "t") {
        P.repAllowed = true;
    }
    else {
        P.repAllowed = false;
    }

    cout << "probablePositive (0 lub 1, domyslnie 0): ";
    getline(cin, linia);
    if (linia.empty()) {
        P.probablePositive = false;
    }
    else {
        int wartosc = stoi(linia);
        P.probablePositive = (wartosc != 0);
    }
}

void zapiszInstancjeDoPliku() {
    if (!I.sekwencjaIstnieje) {
        cout << "Brak instancji do zapisania." << endl;
        return;
    }

    cout << "Podaj nazwe pliku instancji: ";
    string nazwa;
    getline(cin, nazwa);

    ofstream plikWyjsciowy(nazwa);
    if (!plikWyjsciowy) {
        cout << "Blad otwarcia pliku do zapisu." << endl;
        return;
    }
    plikWyjsciowy << I.sekwencja << '\n';
    plikWyjsciowy << P.n << ' '
        << P.k << ' '
        << P.delta_k << ' '
        << P.l_neg << ' '
        << P.l_poz << ' '
        << (P.repAllowed ? 1 : 0) << ' '
        << (P.probablePositive ? 1 : 0) << '\n';
    if (!I.spektrumZBledamiIstnieje) {
        generujSpektrumBezBledow();
        zastosujBledyWSpektrum();
    }

    plikWyjsciowy << I.spektrumZBledami.size() << '\n';
    for (const auto& oligo : I.spektrumZBledami) {
        plikWyjsciowy << oligo << '\n';
    }

    cout << "Instancja zostala zapisana do pliku: " << nazwa << endl;
}

char losujInnyNukleotyd(char aktualny) {
    char nowy;
    do {
        nowy = NUKLEOTYDY[rand() % 4];
    } while (nowy == aktualny);
    return nowy;
}

void generujSpektrumBezBledow() {
    I.spektrum.clear();
    I.spektrumIstnieje = false;

    if (!I.sekwencjaIstnieje) {
        cout << "Brak DNA do wygenerowania spektrum." << endl;
        return;
    }

    if (P.k <= 0 || P.n <= 0 || P.k > P.n) {
        cout << "Niepoprawne parametry n/k do generowania spektrum." << endl;
        return;
    }

    int n = I.sekwencja.size();
    int k = P.k;
    int deltaMax = P.delta_k;

    int ostatniaPozycja = n - k;
    if (ostatniaPozycja < 0) {
        cout << "Sekwencja za krotka do generowania oligonukleotydow." << endl;
        return;
    }

    for (int i = 0; i <= ostatniaPozycja; ++i) {
        int dlugoscOligo = k;

        if (deltaMax > 0) {
            bool moznaZmieniacDlugosc = false;
            int last0 = ostatniaPozycja;          
            int last1 = ostatniaPozycja - 1;      
            int last2 = ostatniaPozycja - 2;      

            if (i != 0 && i != last0 && i != last1 && i != last2) {
                moznaZmieniacDlugosc = true;
            }

            if (moznaZmieniacDlugosc) {
                int delta = rand() % (deltaMax + 1);
                if (delta != 0) {
                    int znak = (rand() % 2 == 0) ? -1 : 1;
                    int potencjalnaDlugosc = k + znak * delta;

                    if (potencjalnaDlugosc < 1 || i + potencjalnaDlugosc > n) {
                        potencjalnaDlugosc = k;
                    }
                    dlugoscOligo = potencjalnaDlugosc;
                }
            }
        }

        string oligo = I.sekwencja.substr(i, dlugoscOligo);

        if (!P.repAllowed) {
            if (find(I.spektrum.begin(), I.spektrum.end(), oligo) != I.spektrum.end()) {
                continue;
            }
        }

        I.spektrum.push_back(oligo);
    }

    I.spektrumIstnieje = !I.spektrum.empty();
}

void zastosujBledyWSpektrum() {
    I.spektrumZBledami.clear();
    I.spektrumZBledamiIstnieje = false;

    if (!I.spektrumIstnieje) {
        cout << "Brak bazowego spektrum do wprowadzenia bledow." << endl;
        return;
    }

    I.spektrumZBledami = I.spektrum;

    int m = I.spektrumZBledami.size();
    if (m == 0) {
        return;
    }

    if (P.l_neg > 0) {
        int maxNeg = m * 0.15;
        int ileNeg = P.l_neg;
        if (maxNeg > 0 && ileNeg > maxNeg) {
            ileNeg = maxNeg;
        }
        if (ileNeg > m) {
            ileNeg = m;
        }

        if (ileNeg > 0) {
            // Chcemy zachowac pierwszy oligo (prefiks), zeby algorytm mial poprawny start.
            int indeksPierwszegoOligo = -1;
            if (!I.pierwszyOligo.empty()) {
                for (int i = 0; i < m; ++i) {
                    if (I.spektrumZBledami[i] == I.pierwszyOligo) {
                        indeksPierwszegoOligo = i;
                        break;
                    }
                }
            } else {
                indeksPierwszegoOligo = 0;
            }

            vector<int> indeksy(m);
            for (int i = 0; i < m; ++i) indeksy[i] = i;
            static random_device rd;
            static mt19937 gen(rd());
            shuffle(indeksy.begin(), indeksy.end(), gen);

            vector<bool> doUsuniecia(m, false);
            int usuniete = 0;
            for (int i = 0; i < m && usuniete < ileNeg; ++i) {
                int idx = indeksy[i];
                if (idx == indeksPierwszegoOligo) {
                    continue;
                }
                doUsuniecia[idx] = true;
                ++usuniete;
            }

            vector<string> bezNeg;
            bezNeg.reserve(m - ileNeg);
            for (int i = 0; i < m; ++i) {
                if (!doUsuniecia[i]) {
                    bezNeg.push_back(I.spektrumZBledami[i]);
                }
            }
            I.spektrumZBledami.swap(bezNeg);
            m = I.spektrumZBledami.size();
        }
    }

    if (P.l_poz > 0) {
        int ilePoz = P.l_poz;

        if (!P.probablePositive) {
            for (int i = 0; i < ilePoz; ++i) {
                int dlugoscOligo = P.k;
                if (P.delta_k > 0) {
                    int delta = rand() % (P.delta_k + 1);
                    if (delta != 0) {
                        int znak = (rand() % 2 == 0) ? -1 : 1;
                        int potencjalnaDlugosc = P.k + znak * delta;
                        if (potencjalnaDlugosc < 1) potencjalnaDlugosc = P.k;
                        dlugoscOligo = potencjalnaDlugosc;
                    }
                }

                string oligo;
                for (int j = 0; j < dlugoscOligo; ++j) {
                    oligo.push_back(NUKLEOTYDY[rand() % 4]);
                }
                if (P.repAllowed || find(I.spektrumZBledami.begin(), I.spektrumZBledami.end(), oligo) == I.spektrumZBledami.end()) {
                    I.spektrumZBledami.push_back(oligo);
                }
            }
        }
        else {
            int dodane = 0;
            int bazowaLiczba = I.spektrum.size();
            if (bazowaLiczba > 0) {
                while (dodane < ilePoz) {
                    const string& bazowy = I.spektrum[rand() % bazowaLiczba];
                    if (bazowy.empty()) {
                        continue;
                    }

                    string e1 = bazowy;
                    e1.back() = losujInnyNukleotyd(e1.back());
                    if (P.repAllowed || find(I.spektrumZBledami.begin(), I.spektrumZBledami.end(), e1) == I.spektrumZBledami.end()) {
                        I.spektrumZBledami.push_back(e1);
                        ++dodane;
                    }
                    if (dodane >= ilePoz) break;

                    string e2 = bazowy;
                    int srodkowyIndex = e2.size() / 2;
                    e2[srodkowyIndex] = losujInnyNukleotyd(e2[srodkowyIndex]);
                    if (P.repAllowed || find(I.spektrumZBledami.begin(), I.spektrumZBledami.end(), e2) == I.spektrumZBledami.end()) {
                        I.spektrumZBledami.push_back(e2);
                        ++dodane;
                    }
                }
            }
        }
    }

    sort(I.spektrumZBledami.begin(), I.spektrumZBledami.end());
    I.spektrumZBledamiIstnieje = !I.spektrumZBledami.empty();
}

void wygenerujSekwencje() {
    cout << "Czy instancja ma byc standardowa? (T/t = tak, inny klawisz = nie): ";
    char odp;
    cin >> odp;
    cin.ignore(numeric_limits<streamsize>::max(), '\n');

    if (odp == 'T' || odp == 't') {
        P.n = 400;
        P.k = 8;
        P.delta_k = 2;
        P.l_neg = 0;
        P.l_poz = 0;
        P.repAllowed = true;
        P.probablePositive = false;
    }
    else {
        pobierzParametryInstancji();
    }
    I.sekwencja.clear();

    for (int i = 0; i < P.n; ++i) {
        I.sekwencja.push_back(NUKLEOTYDY[rand() % 4]);
    }

    I.sekwencjaIstnieje = true;

    generujSpektrumBezBledow();
    // Pierwszy oligo jest zawsze prefiksem DNA (i==0 ma dlugosc k).
    I.pierwszyOligo = (P.k > 0 && I.sekwencja.size() >= P.k)
        ? I.sekwencja.substr(0, P.k)
        : "";
    zastosujBledyWSpektrum();
    

    cout << "\nWygenerowano DNA o dlugosci " << I.sekwencja.size() << endl;
    string podglad = I.sekwencja.substr(0, min<size_t>(80, I.sekwencja.size()));
    cout << "Podglad (pierwsze " << podglad.size() << " znaki): " << podglad << endl;
    cout << "Parametry: n = " << P.n
        << ", k = " << P.k
        << ", delta_k = " << P.delta_k
        << ", l_neg = " << P.l_neg
        << ", l_poz = " << P.l_poz
        << ", repAllowed = " << (P.repAllowed ? "true" : "false")
        << ", probablePositive = " << (P.probablePositive ? 1 : 0)
        << endl;

    if (I.spektrumIstnieje) {
        cout << "Liczba elementow spektrum (bez bledow): " << I.spektrum.size() << endl;
    }
    if (I.spektrumZBledamiIstnieje) {
        cout << "Liczba elementow spektrum z bledami: " << I.spektrumZBledami.size() << endl;
        cout << "Przyklad pierwszych kilku oligo (z bledami):" << endl;
        size_t ile = min<size_t>(5, I.spektrumZBledami.size());
        for (size_t i = 0; i < ile; ++i) {
            cout << i << ": " << I.spektrumZBledami[i]
                 << " (dlugosc = " << I.spektrumZBledami[i].size() << ")" << endl;
        }
    }

    cout << "Czy zapisac instancje do pliku? (T/t = tak, inny klawisz = nie): ";
    char zap;
    cin >> zap;
    cin.ignore(numeric_limits<streamsize>::max(), '\n');

    if (zap == 'T' || zap == 't') {
        zapiszInstancjeDoPliku();
    }
}

void menuGeneratora() {
    char wybor = '9';
    while (wybor != '0') {
        cout << "--- GENERATOR INSTANCJI ---" << endl;
        cout << "1. Wczytaj z pliku" << endl;
        cout << "2. Wygeneruj sekwencje" << endl;
        cout << "0. Powrot" << endl;
        cout << "Wybor: ";
        cin >> wybor;

        switch (wybor) {
        case '1':
            wczytajZPliku();
            break;
        case '2':
            wygenerujSekwencje();
            break;
        case '0':
            break;
        case 'c':
            system("cls");
            break;
        default:
            cout << endl << endl << "Niepoprawny wybor." << endl<< endl << endl;
            system("pause");
            system("cls");
        }
    }
}
int policzNajdluzszyOverlapOgraniczony(const string& aktualnaSekwencja, const string& kandydat) {
    int dlugoscSekwencji = aktualnaSekwencja.size();
    int dlugoscKandydata = kandydat.size();

    int maxOverlap = min({ dlugoscSekwencji, dlugoscKandydata });
    int minOverlap = 1;

    for (int dlugoscOverlap = maxOverlap; dlugoscOverlap >= minOverlap; --dlugoscOverlap) {
        bool pasuje = true;
        for (int i = 0; i < dlugoscOverlap; ++i) {
            if (aktualnaSekwencja[dlugoscSekwencji - dlugoscOverlap + i] != kandydat[i]) {
                pasuje = false;
                break;
            }
        }
        if (pasuje) {
            return dlugoscOverlap;
        }
    }
    return 0;
}

string zlozSekwencjeNaiwnieZachlannie(const Instancja& instancja) {
    if (instancja.spektrumZBledami.empty()) {
        return "";
    }

    int docelowaDlugosc = P.n;

    string aktualnaSekwencja = instancja.pierwszyOligo;

    vector<bool> czyOligoUzyty(instancja.spektrumZBledami.size(), false);
	
    if (!instancja.pierwszyOligo.empty()) {
        for (int i = 0; i < instancja.spektrumZBledami.size(); ++i) {
            if (instancja.spektrumZBledami[i] == instancja.pierwszyOligo) {
                czyOligoUzyty[i] = true;
                break;
            }
        }
    } else {
        czyOligoUzyty[0] = true;
    }

    int limitIteracji = docelowaDlugosc * 2;
    int licznikIteracji = 0;

    while (aktualnaSekwencja.size() < docelowaDlugosc && licznikIteracji < limitIteracji) {
        int najlepszyOverlap = 0;
        int indeksNajlepszego = -1;

        for (int i = 0; i < instancja.spektrumZBledami.size(); ++i) {
            const string& kandydat = instancja.spektrumZBledami[i];
            int overlap = policzNajdluzszyOverlapOgraniczony(
                aktualnaSekwencja,
                kandydat
            );
            if (overlap >= kandydat.size()) {
                continue;
            }
            if (overlap > najlepszyOverlap) {
                najlepszyOverlap = overlap;
                indeksNajlepszego = i;
            }
        }

        if (najlepszyOverlap == 0 || indeksNajlepszego < 0) {
            int indeksNajdluzszego = -1;
            int najlepszaDlugosc = 0;
            for (int i = 0; i < instancja.spektrumZBledami.size(); ++i) {
                if (!P.repAllowed && czyOligoUzyty[i]) {
                    continue;
                }
                int dl = instancja.spektrumZBledami[i].size();
                if (dl > najlepszaDlugosc) {
                    najlepszaDlugosc = dl;
                    indeksNajdluzszego = i;
                }
            }

            if (indeksNajdluzszego < 0 || najlepszaDlugosc <= 0) {
                break;
            }

            aktualnaSekwencja += instancja.spektrumZBledami[indeksNajdluzszego];
            if (!P.repAllowed) {
                czyOligoUzyty[indeksNajdluzszego] = true;
            }
            ++licznikIteracji;
            continue;
        }

        aktualnaSekwencja += instancja.spektrumZBledami[indeksNajlepszego].substr(najlepszyOverlap);
        if (!P.repAllowed) {
            czyOligoUzyty[indeksNajlepszego] = true;
        }

        ++licznikIteracji;
    }

    if (aktualnaSekwencja.size() > docelowaDlugosc) {
        return aktualnaSekwencja.substr(0, docelowaDlugosc);
    }
    return aktualnaSekwencja;
}

int policzOdlegloscLevenshteina(const string& sekwencjaPrawdziwa, const string& sekwencjaOdtworzona) {
    int dlugoscPrawdziwej = sekwencjaPrawdziwa.size();
    int dlugoscOdtworzonej = sekwencjaOdtworzona.size();
    if (dlugoscPrawdziwej == 0) return dlugoscOdtworzonej;
    if (dlugoscOdtworzonej == 0) return dlugoscPrawdziwej;

    vector<int> poprzedniWiersz(dlugoscOdtworzonej + 1, 0);
    vector<int> aktualnyWiersz(dlugoscOdtworzonej + 1, 0);
    for (int j = 0; j <= dlugoscOdtworzonej; ++j) {
        poprzedniWiersz[j] = j;
    }

    for (int i = 1; i <= dlugoscPrawdziwej; ++i) {
        aktualnyWiersz[0] = i;
        for (int j = 1; j <= dlugoscOdtworzonej; ++j) {
            int kosztPodmiany = (sekwencjaPrawdziwa[i - 1] == sekwencjaOdtworzona[j - 1]) ? 0 : 1;
            int kosztUsuniecia = poprzedniWiersz[j] + 1;
            int kosztWstawienia = aktualnyWiersz[j - 1] + 1;
            int kosztZmiany = poprzedniWiersz[j - 1] + kosztPodmiany;
            aktualnyWiersz[j] = min({ kosztUsuniecia, kosztWstawienia, kosztZmiany });
        }
        swap(poprzedniWiersz, aktualnyWiersz);
    }

    return poprzedniWiersz[dlugoscOdtworzonej];
}

pair<int, double> policzPokrycieSpektrum(const vector<string>& spektrumZBledami, const string& sekwencjaOdtworzona) {
    int znalezione = 0;
    int calkowitaLiczba = spektrumZBledami.size();
    for (int i = 0; i < calkowitaLiczba; i++) {
        const string& oligo = spektrumZBledami[i];
        if (!oligo.empty() && sekwencjaOdtworzona.find(oligo) != string::npos) {
            ++znalezione;
        }
    }
    double procent = 0.0;
    if (calkowitaLiczba > 0) {
        procent = (static_cast<double>(znalezione) / calkowitaLiczba) * 100.0;
    }
    return { znalezione, procent };
}

void algorytmNaiwny() {
    if (!I.spektrumZBledamiIstnieje) {
        cout << "Brak spektrum do rekonstrukcji." << endl;
        return;
    }

    string sekwencjaOdtworzona = zlozSekwencjeNaiwnieZachlannie(I);

    cout << "Odtworzona sekwencja: " << endl;
    if (sekwencjaOdtworzona.empty()) {
        cout << "Brak wyniku skladania." << endl;
    } else {
        cout << sekwencjaOdtworzona << endl;
    }

    int odlegloscLevenshteina = 0;
    if (I.sekwencjaIstnieje) {
        odlegloscLevenshteina = policzOdlegloscLevenshteina(I.sekwencja, sekwencjaOdtworzona);
        cout << "Odleglosc Levenshteina: " << odlegloscLevenshteina << endl;
    }

    pair<int, double> wynikPokrycia = policzPokrycieSpektrum(I.spektrumZBledami, sekwencjaOdtworzona);
    cout << "Znalezione elementy spektrum: " << wynikPokrycia.first << "/" << I.spektrumZBledami.size()
    << " (" << wynikPokrycia.second << "% )" << endl;
    pair<int, double> wynikPokryciaOryginalne = policzPokrycieSpektrum(I.spektrum, sekwencjaOdtworzona);
    cout << "Znalezione elementy spektrum: " << wynikPokryciaOryginalne.first << "/" << I.spektrum.size()
         << " (" << wynikPokryciaOryginalne.second << "% )" << endl;

    cout << "Czy zapisac wynik do pliku? (T/N): ";
    char decyzja;
    cin >> decyzja;
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    if (decyzja == 'T' || decyzja == 't') {
        cout << "Podaj nazwe pliku wynikowego: ";
        string nazwaPliku;
        getline(cin, nazwaPliku);
        ofstream plikWyjsciowy(nazwaPliku);
        if (!plikWyjsciowy) {
            cout << "Blad otwarcia pliku do zapisu." << endl;
            return;
        }

        plikWyjsciowy << "Oryginalne DNA:\n" << I.sekwencja << "\n\n";
        plikWyjsciowy << "Odtworzone DNA:\n" << sekwencjaOdtworzona << "\n\n";
        plikWyjsciowy << "Odleglosc Levenshteina: " << odlegloscLevenshteina << "\n";
        plikWyjsciowy << "Pokrycie spektrum: " << wynikPokrycia.first << "/" << I.spektrumZBledami.size()
                      << " (" << wynikPokrycia.second << "% )\n\n";
        plikWyjsciowy << "Parametry instancji:\n";
        plikWyjsciowy << "n = " << P.n << "\n";
        plikWyjsciowy << "k = " << P.k << "\n";
        plikWyjsciowy << "delta_k = " << P.delta_k << "\n";
        plikWyjsciowy << "l_neg = " << P.l_neg << "\n";
        plikWyjsciowy << "l_poz = " << P.l_poz << "\n";
        plikWyjsciowy << "repAllowed = " << (P.repAllowed ? 1 : 0) << "\n";
        plikWyjsciowy << "probablePositive = " << (P.probablePositive ? 1 : 0) << "\n";

        cout << "Wynik zapisany do pliku: " << nazwaPliku << endl;
    }
}

WynikSciezki zbudujRozwiazanieMrowki(
    const vector<string>& wierzcholki,
    const vector<vector<KrawedzGrafu>>& graf,
    vector<vector<double>>& feromony,
    int indeksStartowy,
    int docelowaDlugosc,
    const ParametryMrowkowe& parametry,
    mt19937& generatorLosowy,
    bool uzywajDijkstra
) {
    WynikSciezki wynik;
    wynik.sciezkaWierzcholkow.clear();
    wynik.uzyteKrawedzie.clear();

    int liczbaWierzcholkow = wierzcholki.size();
    vector<int> licznikOdwiedzin(liczbaWierzcholkow, 0);
    vector<char> czyOdwiedzony(liczbaWierzcholkow, 0);

    int obecny = indeksStartowy;
    wynik.sciezkaWierzcholkow.push_back(obecny);
    licznikOdwiedzin[obecny]++;
    czyOdwiedzony[obecny] = 1;
    wynik.sekwencja = wierzcholki[obecny];

    int limitKrokow = docelowaDlugosc * 2;
    int krok = 0;

    while (wynik.sekwencja.size() < docelowaDlugosc && krok < limitKrokow) {
        // Opcjonalny skok Dijkstry co pewna liczbe krokow.
        if (uzywajDijkstra && parametry.coIleKrokowDijkstra > 0 && krok > 0 && (krok % parametry.coIleKrokowDijkstra == 0)) {
            vector<int> nieOdwiedzone;
            nieOdwiedzone.reserve(liczbaWierzcholkow);
            for (int v = 0; v < liczbaWierzcholkow; ++v) {
                if (!czyOdwiedzony[v]) {
                    nieOdwiedzone.push_back(v);
                }
            }
            if (!nieOdwiedzone.empty()) {
                uniform_int_distribution<int> losujCel(0, nieOdwiedzone.size() - 1);
                int cel = nieOdwiedzone[losujCel(generatorLosowy)];
                vector<int> sciezka = dijkstraSciezka(graf, obecny, cel, P.k, P.delta_k);
                if (sciezka.size() >= 2) {
                    // Doklejamy sciezke Dijkstry.
                    for (size_t idx = 1; idx < sciezka.size(); ++idx) {
                        int nastepny = sciezka[idx];
                        int overlap = policzNajdluzszyOverlapDowolny(wynik.sekwencja, wierzcholki[nastepny]);
                        if (overlap >= wierzcholki[nastepny].size()) {
                            continue;
                        }
                        int kategoria = policzKategorieKrawedzi(overlap, P.k, P.delta_k);
                        if (overlap == 0) {
                            wynik.liczbaKrawedziKat3++;
                            wynik.sekwencja += wierzcholki[nastepny];
                        } else {
                            if (kategoria == 1) wynik.liczbaKrawedziKat1++;
                            else if (kategoria == 2) wynik.liczbaKrawedziKat2++;
                            else wynik.liczbaKrawedziKat3++;
                            wynik.sekwencja += wierzcholki[nastepny].substr(overlap);
                        }
                        wynik.uzyteKrawedzie.push_back({ obecny, nastepny });
                        obecny = nastepny;
                        wynik.sciezkaWierzcholkow.push_back(obecny);
                        licznikOdwiedzin[obecny]++;
                        czyOdwiedzony[obecny] = 1;
                        if (wynik.sekwencja.size() >= docelowaDlugosc) {
                            break;
                        }
                    }
                    ++krok;
                    continue;
                }
            }
        }

        // Standardowy krok ACO: wybierz nastepny wierzcholek z rozkladu.
        vector<pair<int, const KrawedzGrafu*>> kandydaci;
        kandydaci.reserve(graf[obecny].size());
        for (const auto& e : graf[obecny]) {
            kandydaci.push_back({ e.doWierzcholka, &e });
        }

        if (!kandydaci.empty()) {
            // Ogranicz liczbe kandydatow (dla wydajnosci): sortujemy po overlap desc.
            sort(kandydaci.begin(), kandydaci.end(), [&](const auto& a, const auto& b) {
                return a.second->overlap > b.second->overlap;
            });
            if (parametry.limitKandydatow > 0 && kandydaci.size() > parametry.limitKandydatow) {
                kandydaci.resize(parametry.limitKandydatow);
            }

            // Jezeli istnieje rozwidlenie z oligo prawie-identycznymi (1-3 roznice),
            // to warto je eksplorowac (typowy efekt bledow dodatnich).
            vector<char> maPodobnego(kandydaci.size(), 0);
            if (kandydaci.size() >= 2) {
                for (size_t i = 0; i < kandydaci.size(); ++i) {
                    const string& si = wierzcholki[kandydaci[i].first];
                    for (size_t j = i + 1; j < kandydaci.size(); ++j) {
                        const string& sj = wierzcholki[kandydaci[j].first];
                        int d = policzOdlegloscHammingaDoLimitu(si, sj, 3);
                        if (d >= 1 && d <= 3) {
                            maPodobnego[i] = 1;
                            maPodobnego[j] = 1;
                        }
                    }
                }
            }

            vector<double> wagi;
            wagi.reserve(kandydaci.size());
            double suma = 0.0;
            for (size_t idx = 0; idx < kandydaci.size(); ++idx) {
                int nastepny = kandydaci[idx].first;
                const KrawedzGrafu& e = *kandydaci[idx].second;

                double tau = feromony[obecny][nastepny];
                double eta = policzHeurystyke(e.overlap, e.kategoriaWagi, czyOdwiedzony[nastepny] != 0);
                if (maPodobnego[idx]) {
                    eta *= 1.35;
                }
                double wartosc = pow(tau, parametry.alfa) * pow(eta, parametry.beta);
                if (wartosc < 1e-12) wartosc = 1e-12;
                wagi.push_back(wartosc);
                suma += wartosc;
            }

            uniform_real_distribution<double> losuj(0.0, suma);
            double prog = losuj(generatorLosowy);
            size_t wybranyIndex = 0;
            double akumulator = 0.0;
            for (size_t i = 0; i < wagi.size(); ++i) {
                akumulator += wagi[i];
                if (akumulator >= prog) {
                    wybranyIndex = i;
                    break;
                }
            }

            int nastepny = kandydaci[wybranyIndex].first;
            const KrawedzGrafu& e = *kandydaci[wybranyIndex].second;
            int overlap = e.overlap;

            if (overlap <= 0) {
                wynik.liczbaKrawedziKat3++;
                wynik.sekwencja += wierzcholki[nastepny];
            } else {
                if (e.kategoriaWagi == 1) wynik.liczbaKrawedziKat1++;
                else if (e.kategoriaWagi == 2) wynik.liczbaKrawedziKat2++;
                else wynik.liczbaKrawedziKat3++;
                wynik.sekwencja += wierzcholki[nastepny].substr(overlap);
            }

            wynik.uzyteKrawedzie.push_back({ obecny, nastepny });
            obecny = nastepny;
            wynik.sciezkaWierzcholkow.push_back(obecny);
            licznikOdwiedzin[obecny]++;
            czyOdwiedzony[obecny] = 1;
        } else {
            // Brak krawedzi z overlap>0: start nowego fragmentu (doklejamy jakis jeszcze nieodwiedzony oligo).
            int wybrany = -1;
            int najlepszaDlugosc = 0;
            for (int v = 0; v < liczbaWierzcholkow; ++v) {
                if (czyOdwiedzony[v]) continue;
                int dl = wierzcholki[v].size();
                if (dl > najlepszaDlugosc) {
                    najlepszaDlugosc = dl;
                    wybrany = v;
                }
            }
            if (wybrany == -1) {
                // Wszystko odwiedzone, bierzemy cokolwiek.
                uniform_int_distribution<int> losujW(0, liczbaWierzcholkow - 1);
                wybrany = losujW(generatorLosowy);
            }

            wynik.liczbaKrawedziKat3++;
            wynik.sekwencja += wierzcholki[wybrany];
            wynik.uzyteKrawedzie.push_back({ obecny, wybrany });
            obecny = wybrany;
            wynik.sciezkaWierzcholkow.push_back(obecny);
            licznikOdwiedzin[obecny]++;
            czyOdwiedzony[obecny] = 1;
        }

        ++krok;
    }

    // Odetnij do docelowej dlugosci.
    if (wynik.sekwencja.size() > docelowaDlugosc) {
        wynik.sekwencja = wynik.sekwencja.substr(0, docelowaDlugosc);
    }

    int unikalne = 0;
    int powtorki = 0;
    for (int v = 0; v < liczbaWierzcholkow; ++v) {
        if (licznikOdwiedzin[v] > 0) {
            ++unikalne;
            if (licznikOdwiedzin[v] > 1) {
                powtorki += (licznikOdwiedzin[v] - 1);
            }
        }
    }
    wynik.unikalneWierzcholki = unikalne;
    wynik.powtorzeniaWierzcholkow = powtorki;
    return wynik;
}

WynikSciezki zbudujRozwiazanieNaiwneDokladnie(
    const Instancja& instancja,
    const vector<string>& wierzcholki,
    int docelowaDlugosc
) {
    WynikSciezki wynik;
    wynik.sciezkaWierzcholkow.clear();
    wynik.uzyteKrawedzie.clear();

    if (wierzcholki.empty() || docelowaDlugosc <= 0) {
        wynik.sekwencja = "";
        return wynik;
    }

    string aktualnaSekwencja;
    int obecnyW = 0;
    if (!instancja.pierwszyOligo.empty()) {
        int idx = znajdzIndeksWierzcholka(wierzcholki, instancja.pierwszyOligo);
        if (idx >= 0) {
            obecnyW = idx;
            aktualnaSekwencja = instancja.pierwszyOligo;
        }
    }
    if (aktualnaSekwencja.empty()) {
        obecnyW = 0;
        aktualnaSekwencja = wierzcholki[0];
    }

    wynik.sciezkaWierzcholkow.push_back(obecnyW);
    wynik.sekwencja = aktualnaSekwencja;

    vector<bool> czyOligoUzyty(wierzcholki.size(), false);
    czyOligoUzyty[obecnyW] = true;

    int limitIteracji = docelowaDlugosc * 2;
    int licznikIteracji = 0;

    while (aktualnaSekwencja.size() < docelowaDlugosc && licznikIteracji < limitIteracji) {
        int najlepszyOverlap = 0;
        int indeksNajlepszego = -1;

        for (size_t i = 0; i < wierzcholki.size(); ++i) {
            if (!P.repAllowed && czyOligoUzyty[i]) {
                continue;
            }

            const string& kandydat = wierzcholki[i];
            int overlap = policzNajdluzszyOverlapOgraniczony(aktualnaSekwencja, kandydat);

            if (overlap >= kandydat.size()) {
                continue;
            }

            if (overlap > najlepszyOverlap) {
                najlepszyOverlap = overlap;
                indeksNajlepszego = i;
            }
        }

        if (najlepszyOverlap == 0 || indeksNajlepszego < 0) {
            int indeksNajdluzszego = -1;
            int najlepszaDlugosc = 0;
            for (size_t i = 0; i < wierzcholki.size(); ++i) {
                if (!P.repAllowed && czyOligoUzyty[i]) {
                    continue;
                }
                int dl = wierzcholki[i].size();
                if (dl > najlepszaDlugosc) {
                    najlepszaDlugosc = dl;
                    indeksNajdluzszego = i;
                }
            }

            if (indeksNajdluzszego < 0 || najlepszaDlugosc <= 0) {
                break;
            }

            const string& nastepnyOligo = wierzcholki[indeksNajdluzszego];
            aktualnaSekwencja += nastepnyOligo;
            if (!P.repAllowed) {
                czyOligoUzyty[indeksNajdluzszego] = true;
            }

            wynik.liczbaKrawedziKat3++;
            wynik.uzyteKrawedzie.push_back({ obecnyW, indeksNajdluzszego });
            obecnyW = indeksNajdluzszego;
            wynik.sciezkaWierzcholkow.push_back(obecnyW);
            ++licznikIteracji;
            continue;
        }

        const string& nastepnyOligo = wierzcholki[indeksNajlepszego];
        aktualnaSekwencja += nastepnyOligo.substr(najlepszyOverlap);
        if (!P.repAllowed) {
            czyOligoUzyty[indeksNajlepszego] = true;
        }

        int kategoria = policzKategorieKrawedzi(najlepszyOverlap, P.k, P.delta_k);
        if (kategoria == 1) wynik.liczbaKrawedziKat1++;
        else if (kategoria == 2) wynik.liczbaKrawedziKat2++;
        else wynik.liczbaKrawedziKat3++;
        wynik.uzyteKrawedzie.push_back({ obecnyW, indeksNajlepszego });
        obecnyW = indeksNajlepszego;
        wynik.sciezkaWierzcholkow.push_back(obecnyW);

        ++licznikIteracji;
    }

    if (aktualnaSekwencja.size() > docelowaDlugosc) {
        aktualnaSekwencja = aktualnaSekwencja.substr(0, docelowaDlugosc);
    }
    wynik.sekwencja = aktualnaSekwencja;

    if (!wierzcholki.empty()) {
        vector<int> licznikOdwiedzin(wierzcholki.size(), 0);
        for (int v : wynik.sciezkaWierzcholkow) {
            if (v >= 0 && v < licznikOdwiedzin.size()) {
                licznikOdwiedzin[v]++;
            }
        }
        int unikalne = 0;
        int powtorki = 0;
        for (int c : licznikOdwiedzin) {
            if (c > 0) {
                ++unikalne;
                if (c > 1) powtorki += (c - 1);
            }
        }
        wynik.unikalneWierzcholki = unikalne;
        wynik.powtorzeniaWierzcholkow = powtorki;
    }

    return wynik;
}

WynikSciezki zbudujRozwiazanieNaiwneDlaAco(
    const vector<string>& wierzcholki,
    const vector<vector<KrawedzGrafu>>& graf,
    int indeksStartowy,
    int docelowaDlugosc,
    int k,
    int delta_k
) {
    WynikSciezki wynik;
    wynik.sciezkaWierzcholkow.clear();
    wynik.uzyteKrawedzie.clear();

    int liczbaWierzcholkow = wierzcholki.size();
    vector<int> licznikOdwiedzin(liczbaWierzcholkow, 0);
    vector<char> czyOdwiedzony(liczbaWierzcholkow, 0);

    int obecny = indeksStartowy;
    wynik.sciezkaWierzcholkow.push_back(obecny);
    licznikOdwiedzin[obecny]++;
    czyOdwiedzony[obecny] = 1;
    wynik.sekwencja = wierzcholki[obecny];

    int limitKrokow = docelowaDlugosc * 2;
    int krok = 0;
    int maxOverlap = k + delta_k;
    int minOverlap = std::max(1, k - delta_k);

    while (wynik.sekwencja.size() < docelowaDlugosc && krok < limitKrokow) {
        int wybrany = -1;
        int najlepszyOverlap = 0;
        int najlepszaKategoria = 99;
        bool najlepszyNieodwiedzony = false;

        // Dwie fazy jak w naiwnym:
        // 1) overlap w zakresie (k ± delta_k)
        // 2) jesli brak - luzujemy i bierzemy dowolny overlap >= 1
        for (int proba = 0; proba < 2 && wybrany == -1; proba++) {
            for (const auto& e : graf[obecny]) {
                int nastepny = e.doWierzcholka;
                if (e.overlap <= 0) continue;
                if (e.overlap >= wierzcholki[nastepny].size()) continue;

                if (proba == 0) {
                    if (e.overlap < minOverlap || e.overlap > maxOverlap) continue;
                }

                bool nieodwiedzony = (licznikOdwiedzin[nastepny] == 0);

                bool lepszy = false;
                if (e.overlap > najlepszyOverlap) lepszy = true;
                else if (e.overlap == najlepszyOverlap) {
                    if (nieodwiedzony && !najlepszyNieodwiedzony) lepszy = true;
                    else if (nieodwiedzony == najlepszyNieodwiedzony && e.kategoriaWagi < najlepszaKategoria) lepszy = true;
                }

                if (lepszy) {
                    najlepszyOverlap = e.overlap;
                    najlepszaKategoria = e.kategoriaWagi;
                    najlepszyNieodwiedzony = nieodwiedzony;
                    wybrany = nastepny;
                }
            }
        }

        if (wybrany != -1 && najlepszyOverlap > 0) {
            if (najlepszaKategoria == 1) wynik.liczbaKrawedziKat1++;
            else if (najlepszaKategoria == 2) wynik.liczbaKrawedziKat2++;
            else wynik.liczbaKrawedziKat3++;
            wynik.sekwencja += wierzcholki[wybrany].substr(najlepszyOverlap);
            wynik.uzyteKrawedzie.push_back({ obecny, wybrany });
            obecny = wybrany;
            wynik.sciezkaWierzcholkow.push_back(obecny);
            licznikOdwiedzin[obecny]++;
            czyOdwiedzony[obecny] = 1;
        } else {
            // Brak ruchu: tak jak w naiwnym - doklejamy nowy fragment (najdluzszy nieodwiedzony).
            int indeksNajdluzszego = -1;
            int najlepszaDlugosc = 0;
            for (int v = 0; v < liczbaWierzcholkow; ++v) {
                if (czyOdwiedzony[v]) continue;
                int dl = wierzcholki[v].size();
                if (dl > najlepszaDlugosc) {
                    najlepszaDlugosc = dl;
                    indeksNajdluzszego = v;
                }
            }
            if (indeksNajdluzszego < 0) {
                break;
            }
            wynik.liczbaKrawedziKat3++;
            wynik.sekwencja += wierzcholki[indeksNajdluzszego];
            wynik.uzyteKrawedzie.push_back({ obecny, indeksNajdluzszego });
            obecny = indeksNajdluzszego;
            wynik.sciezkaWierzcholkow.push_back(obecny);
            licznikOdwiedzin[obecny]++;
            czyOdwiedzony[obecny] = 1;
        }

        ++krok;
    }

    if (wynik.sekwencja.size() > docelowaDlugosc) {
        wynik.sekwencja = wynik.sekwencja.substr(0, docelowaDlugosc);
    }

    int unikalne = 0;
    int powtorki = 0;
    for (int v = 0; v < liczbaWierzcholkow; v++) {
        if (licznikOdwiedzin[v] > 0) {
            ++unikalne;
            if (licznikOdwiedzin[v] > 1) {
                powtorki += (licznikOdwiedzin[v] - 1);
            }
        }
    }
    wynik.unikalneWierzcholki = unikalne;
    wynik.powtorzeniaWierzcholkow = powtorki;
    return wynik;
}

void metaheurystyka() {
    ParametryMrowkowe parametry;
    cout << "--- METAHEURYSTYKA (Algorytm mrowkowy) ---" << endl;
    cout << "Wcisnij ENTER aby uzyc domyslnych wartosci." << endl;
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    parametry.liczbaIteracji = wczytajIntZDomyslna("Liczba iteracji (domyslnie 300): ", parametry.liczbaIteracji);
    parametry.liczbaMrowek = wczytajIntZDomyslna("Liczba mrowek (domyslnie 150): ", parametry.liczbaMrowek);
    parametry.limitKandydatow = wczytajIntZDomyslna("Limit kandydatow na krok (domyslnie 120): ", parametry.limitKandydatow);
    parametry.procentMrowekDijkstra = wczytajIntZDomyslna("Procent mrowek z Dijkstra (0-100, domyslnie 25): ", parametry.procentMrowekDijkstra);
    parametry.coIleKrokowDijkstra = wczytajIntZDomyslna("Co ile krokow uzyc Dijkstra (domyslnie 50): ", parametry.coIleKrokowDijkstra);

    if (parametry.liczbaIteracji < 1) parametry.liczbaIteracji = 1;
    if (parametry.liczbaMrowek < 1) parametry.liczbaMrowek = 1;
    if (parametry.limitKandydatow < 0) parametry.limitKandydatow = 5;
    if (parametry.procentMrowekDijkstra < 0) parametry.procentMrowekDijkstra = 0;
    if (parametry.procentMrowekDijkstra > 100) parametry.procentMrowekDijkstra = 100;
    if (parametry.coIleKrokowDijkstra < 1) parametry.coIleKrokowDijkstra = 50;

    vector<string> wierzcholki = I.spektrumZBledami;

    wierzcholki.erase(unique(wierzcholki.begin(), wierzcholki.end()), wierzcholki.end());
    if (wierzcholki.empty()) {
        cout << "Spektrum jest puste." << endl;
        return;
    }

    int indeksStartowy = 0;
    int index = znajdzIndeksWierzcholka(wierzcholki, I.pierwszyOligo);
    if (index >= 0) indeksStartowy = index;

    int liczbaWierzcholkow = wierzcholki.size();
    vector<vector<KrawedzGrafu>> graf(liczbaWierzcholkow);
    for (int i = 0; i < liczbaWierzcholkow; i++) {
        for (int j = 0; j < liczbaWierzcholkow; j++) {
            if (i == j) continue;
            int overlap = policzNajdluzszyOverlapDowolny(wierzcholki[i], wierzcholki[j]);
            if (overlap <= 0) continue;
            if (overlap >= wierzcholki[j].size()) continue; // nic nie wniesie
            KrawedzGrafu e;
            e.doWierzcholka = j;
            e.overlap = overlap;
            e.kategoriaWagi = policzKategorieKrawedzi(overlap, P.k, P.delta_k);
            graf[i].push_back(e);
        }
        sort(graf[i].begin(), graf[i].end(), [](const KrawedzGrafu& a, const KrawedzGrafu& b) {
            if (a.kategoriaWagi != b.kategoriaWagi) return a.kategoriaWagi < b.kategoriaWagi;
            return a.overlap > b.overlap;
        });
    }

    // Feromony jako macierz (dla prostoty). Inicjalizacja mala dodatnia.
    vector<vector<double>> feromony(liczbaWierzcholkow, vector<double>(liczbaWierzcholkow, 0.01));

    static random_device rd;
    mt19937 generatorLosowy(rd());
    uniform_int_distribution<int> losujProcent(0, 99);

    WynikSciezki najlepszyGlobalny;
    double najlepszaOcenaGlobalna = -1e18;

    // Pierwsza mrowka ma byc dokladnie algorytmem naiwnym (ta sama logika wyboru oligo).
    // Liczymy raz i potem tylko kopiujemy do iteracji.
    WynikSciezki bazowyNaiwny = zbudujRozwiazanieNaiwneDokladnie(I, wierzcholki, P.n);
    najlepszyGlobalny = bazowyNaiwny;
    najlepszaOcenaGlobalna = policzOceneSciezki(bazowyNaiwny);
    auto najlepszePokrycieGlobalne = policzPokrycieSpektrum(I.spektrumZBledami, najlepszyGlobalny.sekwencja);

    for (int iter = 0; iter < parametry.liczbaIteracji; ++iter) {
        WynikSciezki najlepszyIteracji;
        double najlepszaOcenaIteracji = -1e18;

        for (int mrowka = 0; mrowka < parametry.liczbaMrowek; ++mrowka) {
            WynikSciezki wynik;
            if (mrowka == 0) {
                // Pierwsza mrowka: dokladnie algorytm naiwny.
                wynik = bazowyNaiwny;
            } else {
                bool uzywajDijkstra = (parametry.procentMrowekDijkstra > 0) && (losujProcent(generatorLosowy) < parametry.procentMrowekDijkstra);
                wynik = zbudujRozwiazanieMrowki(
                    wierzcholki,
                    graf,
                    feromony,
                    indeksStartowy,
                    P.n,
                    parametry,
                    generatorLosowy,
                    uzywajDijkstra
                );
            }

            double ocena = policzOceneSciezki(wynik);
            auto pokrycie = policzPokrycieSpektrum(I.spektrumZBledami, wynik.sekwencja);
            if (ocena > najlepszaOcenaIteracji) {
                najlepszaOcenaIteracji = ocena;
                najlepszyIteracji = wynik;
            }
            // Wynik koncowy wybieramy glownie po pokryciu spektrum (bez uzycia oryginalnego DNA),
            // a ocena jest tylko tie-breakerem.
            if (pokrycie.first > najlepszePokrycieGlobalne.first ||
                (pokrycie.first == najlepszePokrycieGlobalne.first && ocena > najlepszaOcenaGlobalna)) {
                najlepszePokrycieGlobalne = pokrycie;
                najlepszaOcenaGlobalna = ocena;
                najlepszyGlobalny = wynik;
            }
        }

        // Parowanie.
        double wspolczynnikParowania = 1.0 - parametry.parowanie;
        for (int i = 0; i < liczbaWierzcholkow; ++i) {
            for (int j = 0; j < liczbaWierzcholkow; ++j) {
                feromony[i][j] *= wspolczynnikParowania;
                if (feromony[i][j] < 1e-6) feromony[i][j] = 1e-6;
            }
        }

        // Depozyt feromonu na najlepszej sciezce iteracji.
        double ocenaDoWzmocnienia = std::max(1.0, najlepszaOcenaIteracji);
        double depozyt = parametry.q * (ocenaDoWzmocnienia / (najlepszyIteracji.uzyteKrawedzie.size() + 1.0));
        for (const auto& e : najlepszyIteracji.uzyteKrawedzie) {
            feromony[e.first][e.second] += depozyt;
        }

        if ((iter + 1) % 20 == 0 || iter == 0 || iter + 1 == parametry.liczbaIteracji) {
            cout << "Iteracja " << (iter + 1) << "/" << parametry.liczbaIteracji
                 << " | ocena(best)=" << najlepszaOcenaGlobalna
                 << " | kat1=" << najlepszyGlobalny.liczbaKrawedziKat1
                 << " kat2=" << najlepszyGlobalny.liczbaKrawedziKat2
                 << " kat3=" << najlepszyGlobalny.liczbaKrawedziKat3
                 << " | unikalneW=" << najlepszyGlobalny.unikalneWierzcholki
                 << " | pokrycie=" << najlepszePokrycieGlobalne.first << "/" << I.spektrumZBledami.size()
                 << endl;
        }
    }

    string sekwencjaOdtworzona = najlepszyGlobalny.sekwencja;
    cout << "Odtworzona sekwencja (dlugosc = " << sekwencjaOdtworzona.size() << "):" << endl;
    cout << sekwencjaOdtworzona << endl;

    int odlegloscLevenshteina = 0;
    if (I.sekwencjaIstnieje) {
        odlegloscLevenshteina = policzOdlegloscLevenshteina(I.sekwencja, sekwencjaOdtworzona);
        cout << "Odleglosc Levenshteina: " << odlegloscLevenshteina << endl;
    }

    auto wynikPokrycia = policzPokrycieSpektrum(I.spektrumZBledami, sekwencjaOdtworzona);
    cout << "Znalezione elementy spektrum: " << wynikPokrycia.first << "/" << I.spektrumZBledami.size()
         << " (" << wynikPokrycia.second << "% )" << endl;
    auto wynikPokryciaOryginalne = policzPokrycieSpektrum(I.spektrum, sekwencjaOdtworzona);
    cout << "Znalezione elementy spektrum: " << wynikPokryciaOryginalne.first << "/" << I.spektrum.size()
         << " (" << wynikPokryciaOryginalne.second << "% )" << endl;

    cout << "Czy zapisac wynik do pliku? (T/N): ";
    char decyzja;
    cin >> decyzja;
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    if (decyzja == 'T' || decyzja == 't') {
        cout << "Podaj nazwe pliku wynikowego: ";
        string nazwaPliku;
        getline(cin, nazwaPliku);
        zapiszWynikDoPliku(nazwaPliku, sekwencjaOdtworzona, odlegloscLevenshteina, wynikPokrycia);
        cout << "Wynik zapisany do pliku: " << nazwaPliku << endl;
    }
}

// ----------------- Tryb automatycznych testow (batch) -----------------
// Uruchamianie:
//   ./main --batch test.txt
// Parametry batch sa zgodne z prosba: 100 iteracji, mrówki 10..100 co 10,
// limit kandydatow=2, procentDijkstra=25, coIleDijkstra=50, 3 proby na punkt.
static bool wczytajInstancjeZPlikuBezInterakcji(const string& nazwaPliku) {
    ifstream in(nazwaPliku);
    if (!in) {
        cerr << "Blad otwarcia pliku: " << nazwaPliku << "\n";
        return false;
    }

    I.sekwencja.clear();

    string dnaLinia;
    in >> ws;
    if (!getline(in, dnaLinia) || dnaLinia.empty()) {
        cerr << "Plik pusty lub blad formatu.\n";
        return false;
    }

    I.sekwencja = dnaLinia;
    I.sekwencjaIstnieje = true;

    string paramLinia;
    if (getline(in, paramLinia) && !paramLinia.empty()) {
        istringstream iss(paramLinia);
        int n, k, d, lneg, lpoz, rep, prob;
        if (iss >> n >> k >> d >> lneg >> lpoz >> rep >> prob) {
            P.n = n;
            P.k = k;
            P.delta_k = d;
            P.l_neg = lneg;
            P.l_poz = lpoz;
            P.repAllowed = (rep != 0);
            P.probablePositive = (prob != 0);
        } else {
            ustawParametryDomyslneDlaAktualnejSekwencji();
        }
    } else {
        ustawParametryDomyslneDlaAktualnejSekwencji();
    }

    I.pierwszyOligo = (P.k > 0 && I.sekwencja.size() >= static_cast<size_t>(P.k))
        ? I.sekwencja.substr(0, P.k)
        : "";

    generujSpektrumBezBledow();
    I.spektrumZBledami.clear();
    I.spektrumZBledamiIstnieje = false;

    int liczbaOligoZPliku = 0;
    if (in >> liczbaOligoZPliku && liczbaOligoZPliku > 0) {
        string usuniecieZnaku;
        getline(in, usuniecieZnaku);
        for (int i = 0; i < liczbaOligoZPliku; ++i) {
            string oligo;
            if (!getline(in, oligo)) break;
            if (!oligo.empty()) {
                I.spektrumZBledami.push_back(oligo);
            }
        }
        sort(I.spektrumZBledami.begin(), I.spektrumZBledami.end());
        I.spektrumZBledamiIstnieje = !I.spektrumZBledami.empty();
    }

    if (!I.spektrumZBledamiIstnieje) {
        zastosujBledyWSpektrum();
    }

    return I.spektrumZBledamiIstnieje;
}

static int uruchomBatchZPliku(const string& nazwaPliku) {
    if (!wczytajInstancjeZPlikuBezInterakcji(nazwaPliku)) {
        return 2;
    }

    ParametryMrowkowe parametry;
    parametry.liczbaIteracji = 100;
    parametry.limitKandydatow = 2;
    parametry.procentMrowekDijkstra = 25;
    parametry.coIleKrokowDijkstra = 50;

    vector<string> wierzcholki = I.spektrumZBledami;
    wierzcholki.erase(unique(wierzcholki.begin(), wierzcholki.end()), wierzcholki.end());
    if (wierzcholki.empty()) {
        cerr << "Spektrum jest puste.\n";
        return 3;
    }

    int indeksStartowy = 0;
    int index = znajdzIndeksWierzcholka(wierzcholki, I.pierwszyOligo);
    if (index >= 0) indeksStartowy = index;

    int liczbaWierzcholkow = static_cast<int>(wierzcholki.size());
    vector<vector<KrawedzGrafu>> graf(liczbaWierzcholkow);
    for (int i = 0; i < liczbaWierzcholkow; i++) {
        for (int j = 0; j < liczbaWierzcholkow; j++) {
            if (i == j) continue;
            int overlap = policzNajdluzszyOverlapDowolny(wierzcholki[i], wierzcholki[j]);
            if (overlap <= 0) continue;
            if (overlap >= static_cast<int>(wierzcholki[j].size())) continue;
            KrawedzGrafu e;
            e.doWierzcholka = j;
            e.overlap = overlap;
            e.kategoriaWagi = policzKategorieKrawedzi(overlap, P.k, P.delta_k);
            graf[i].push_back(e);
        }
        sort(graf[i].begin(), graf[i].end(), [](const KrawedzGrafu& a, const KrawedzGrafu& b) {
            if (a.kategoriaWagi != b.kategoriaWagi) return a.kategoriaWagi < b.kategoriaWagi;
            return a.overlap > b.overlap;
        });
    }

    static random_device rd;
    mt19937 generatorLosowy(rd());
    uniform_int_distribution<int> losujProcent(0, 99);

    cout << "mrowki\tsrednie_pokrycie_procent" << endl;
    for (int mrowki = 10; mrowki <= 100; mrowki += 10) {
        parametry.liczbaMrowek = mrowki;

        double sumaPokryc = 0.0;
        for (int proba = 0; proba < 3; ++proba) {
            vector<vector<double>> feromony(liczbaWierzcholkow, vector<double>(liczbaWierzcholkow, 0.01));

            WynikSciezki najlepszyGlobalny;
            double najlepszaOcenaGlobalna = -1e18;
            auto najlepszePokrycieGlobalne = make_pair(0, 0.0);

            for (int iter = 0; iter < parametry.liczbaIteracji; ++iter) {
                WynikSciezki najlepszyIteracji;
                double najlepszaOcenaIteracji = -1e18;

                for (int mrowka = 0; mrowka < parametry.liczbaMrowek; ++mrowka) {
                    bool uzywajDijkstra = (parametry.procentMrowekDijkstra > 0) && (losujProcent(generatorLosowy) < parametry.procentMrowekDijkstra);
                    WynikSciezki wynik = zbudujRozwiazanieMrowki(
                        wierzcholki,
                        graf,
                        feromony,
                        indeksStartowy,
                        P.n,
                        parametry,
                        generatorLosowy,
                        uzywajDijkstra
                    );

                    double ocena = policzOceneSciezki(wynik);
                    auto pokrycie = policzPokrycieSpektrum(I.spektrumZBledami, wynik.sekwencja);

                    if (ocena > najlepszaOcenaIteracji) {
                        najlepszaOcenaIteracji = ocena;
                        najlepszyIteracji = wynik;
                    }

                    if (pokrycie.first > najlepszePokrycieGlobalne.first ||
                        (pokrycie.first == najlepszePokrycieGlobalne.first && ocena > najlepszaOcenaGlobalna)) {
                        najlepszePokrycieGlobalne = pokrycie;
                        najlepszaOcenaGlobalna = ocena;
                        najlepszyGlobalny = wynik;
                    }
                }

                double wspolczynnikParowania = 1.0 - parametry.parowanie;
                for (int i = 0; i < liczbaWierzcholkow; ++i) {
                    for (int j = 0; j < liczbaWierzcholkow; ++j) {
                        feromony[i][j] *= wspolczynnikParowania;
                        if (feromony[i][j] < 1e-6) feromony[i][j] = 1e-6;
                    }
                }

                double ocenaDoWzmocnienia = std::max(1.0, najlepszaOcenaIteracji);
                double depozyt = parametry.q * (ocenaDoWzmocnienia / (najlepszyIteracji.uzyteKrawedzie.size() + 1.0));
                for (const auto& e : najlepszyIteracji.uzyteKrawedzie) {
                    feromony[e.first][e.second] += depozyt;
                }
            }

            sumaPokryc += najlepszePokrycieGlobalne.second;
        }

        double srednie = sumaPokryc / 3.0;
        cout << mrowki << "\t" << srednie << endl;
    }

    return 0;
}

int main(int argc, char** argv) {
    srand(time(nullptr));

    // Tryb batch: ./main --batch test.txt
    if (argc >= 3) {
        string a1 = argv[1];
        if (a1 == "--batch") {
            return uruchomBatchZPliku(argv[2]);
        }
    }

    char wybor = '1';
    while (wybor != '0') {
        wyswietlanieMenu();
        cout << "Wybor: ";
        cin >> wybor;

        if ((wybor == '2' || wybor == '3') && !I.sekwencjaIstnieje) {
            system("cls"); 
            cout << "Najpierw wybierz/wygeneruj sekwencje." << endl;
            cout<< "-------------------------------------" << endl;
            system("pause"); 
            continue;        
        }

        switch (wybor) {
        case '1':
            menuGeneratora();
            break;
        case '2':
            algorytmNaiwny();
            break;
        case '3':
            metaheurystyka();
            break;
        case 'c': 
            system("cls"); 
            break;
        case '0':
            cout << "Koniec programu." << endl;
            break;
        default:
            cout << endl << endl << "Niepoprawny wybor." << endl<< endl << endl;
            system("pause"); 
            system("cls"); 
            break;
        }
    }
    return 0;
}