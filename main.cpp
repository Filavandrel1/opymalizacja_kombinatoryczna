#include <iostream>
#include <fstream>
#include <string>
using namespace std;

string sekwencja;
bool sekwencjaIstnieje = false;


void wyswietlanieMenu() {
    if (sekwencjaIstnieje) {
        cout << "Aktualna sekwencja: " << sekwencja << endl;
        cout << "Dlugosc: " << sekwencja.size() << endl;
    }
    else {
        cout << "Brak wybranej sekwencji." << endl;
    }
    cout << "------MENU----------" << endl;
    cout << "1. Generator instancji" << endl;
    cout << "2. Algorytm naiwny" << endl;
    cout << "3. Metaheurystyka" << endl;
    cout << "0. Wyjscie z programu" << endl;
    cout << "--------------------" << endl;
}

void zapisSekwencjiDoPliku() {
    if (!sekwencjaIstnieje) {
        cout << "Brak sekwencji do zapisania." << endl;
        return;
    }
    cout << "Podaj nazwe pliku do zapisu: ";
    string nazwa;
    cin >> nazwa;
    ofstream out(nazwa);
    if (!out) {
        cout << "Blad otwarcia pliku." << endl;
        return;
    }
    out << sekwencja;
    cout << "Sekwencja zapisana do pliku." << endl;
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
    sekwencja.clear();
    in >> sekwencja;
    if (!sekwencja.empty()) {
        sekwencjaIstnieje = true;
        cout << "Wczytano sekwencje: " << sekwencja << endl;
    }
    else {
        cout << "Plik pusty lub blad formatu." << endl;
    }
}

void wygenerujSekwencje() {}

void menuGeneratora() {
    int wybor = -1;
    while (wybor != 0) {
        cout << "--- GENERATOR INSTANCJI ---" << endl;
        cout << "1. Wczytaj z pliku" << endl;
        cout << "2. Wygeneruj sekwencje" << endl;
        cout << "0. Powrot" << endl;
        cout << "Wybor: ";
        cin >> wybor;

        switch (wybor) {
        case 1:
            wczytajZPliku();
            break;
        case 2:
            wygenerujSekwencje();
            break;
        case 0:
            break;
        default:
            cout << "Niepoprawny wybor." << endl;
        }
    }
}

void algorytmNaiwny() {}
void metaheurystyka() {}

int main() {
    int wybor = -1;
    while (wybor != 0) {
        wyswietlanieMenu();
        cout << "Wybor: ";
        cin >> wybor;

        switch (wybor) {
        case 1:
            menuGeneratora();
            break;
        case 2:
            if (!sekwencjaIstnieje) {
                cout << "Najpierw wybierz/wygeneruj sekwencje." << endl;
            }
            else {
                algorytmNaiwny();
            }
            break;
        case 3:
            if (!sekwencjaIstnieje) {
                cout << "Najpierw wybierz/wygeneruj sekwencje." << endl;
            }
            else {
                metaheurystyka();
            }
            break;
        case 0:
            cout << "Koniec programu." << endl;
            break;
        default:
            cout << "Niepoprawny wybor." << endl;
        }
    }

    return 0;
}