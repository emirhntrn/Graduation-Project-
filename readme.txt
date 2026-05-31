TR

PROJE ADI
Akışkanla Temas Halindeki Takviyeli Plakaların Çok Amaçlı Yapısal Optimizasyonu

HAZIRLAYAN
Emir Han TURAN

DANIŞMAN
Dr. Mustafa Erden YILDIZDAĞ

KURUM
İstanbul Teknik Üniversitesi
Gemi İnşaatı ve Deniz Bilimleri Fakültesi
Gemi İnşaatı ve Gemi Makineleri Mühendisliği Bölümü

PROJE AÇIKLAMASI
Bu proje kapsamında, akışkanla temas halinde bulunan düz takviyeli bir plaka yapısının ıslak modal davranışı ve yapısal kütlesi incelenmiştir. Yapısal plaka ve takviye profilleri SHELL181 elemanlarıyla, akışkan hacmi ise FLUID30 elemanlarıyla modellenmiştir. Sonlu eleman analizleri ANSYS MAPDL / PyMAPDL ortamında yürütülmüş, çok amaçlı optimizasyon işlemi ise Python ortamında pymoo kütüphanesi ile NSGA-II algoritması kullanılarak gerçekleştirilmiştir.

DOSYA İÇERİĞİ

1. Solver.py
   Parametrik sonlu eleman modelini oluşturan ana çözüm dosyasıdır. Bu dosyada geometri, malzeme özellikleri, eleman tipleri, mesh ayarları, sınır koşulları, FSI tanımı ve ıslak modal analiz adımları yer almaktadır. evaluate_design(...) fonksiyonu, verilen tasarım değişkenleri için seçilen ıslak doğal frekans, yapısal hacim ve yapısal kütle değerlerini döndürür.

2. Test.py
   Solver dosyasının doğru çalışıp çalışmadığını kontrol etmek için hazırlanmıştır. Referans bir tasarım için tek analiz çalıştırır ve frekans, hacim, kütle, profil sayısı ve mod bilgilerini ekrana yazdırır. Optimizasyona başlamadan önce bu dosyanın çalıştırılması önerilir.

3. Optimization.py
   NSGA-II optimizasyon kodudur. Plaka kalınlığı, profil kalınlığı, profil yüksekliği ve profil aralığı değişkenleri üzerinden optimizasyon yapar. Amaçlar, seçilen ıslak doğal frekansın artırılması ve yapısal kütlenin azaltılmasıdır. Optimizasyon sonunda Pareto sonuçları, tüm değerlendirmeler ve seçilen optimum tasarım CSV dosyaları olarak kaydedilir.

4. MAPDL macro dosyası
   Bu dosya zorunlu değildir. ANSYS MAPDL ortamında modelin görselleştirilmesi, mesh kalitesinin kontrol edilmesi ve mod şekillerinin görüntülenmesi için kullanılabilir. Bu macro, çözümden ziyade rapor görselleri, fluid/structure mesh görünümü ve deforme mod şekli çıktıları almak amacıyla hazırlanmıştır.

GEREKLİ KÜTÜPHANELER
Python ortamında aşağıdaki kütüphanelerin kurulu olması gerekir:

* numpy
* pandas
* matplotlib
* ansys-mapdl-core
* pymoo

Örnek kurulum komutu:
pip install numpy pandas matplotlib ansys-mapdl-core pymoo

ANSYS GEREKSİNİMİ
Kodların çalışması için bilgisayarda ANSYS MAPDL kurulu olmalıdır. Çalışmada ANSYS 2023 R1 ortamı kullanılmıştır. PyMAPDL, Python üzerinden MAPDL oturumu başlatmak ve analizleri otomatik yürütmek için kullanılmaktadır.

KULLANIM SIRASI

1. Tüm Python dosyalarını aynı klasöre koyun.
2. Solver.py dosyasının import edilebilir durumda olduğundan emin olun.
3. Önce Test.py dosyasını çalıştırarak tek analiz testini yapın.
4. Test başarılı olursa Optimization.py dosyasını çalıştırarak NSGA-II optimizasyonunu başlatın.
5. Analiz tamamlandıktan sonra oluşan çalışma klasöründe CSV sonuç dosyalarını ve grafik çıktılarını kontrol edin.
6. MAPDL macro dosyası kullanılacaksa ANSYS MAPDL içinde çalıştırılarak mesh ve mod şekli görselleri alınabilir.

NSGA-II AYARLARINI DEĞİŞTİRME
Optimizasyonda farklı sonuçlar almak için Optimization.py dosyasındaki USER SETTINGS bölümündeki NSGA-II parametreleri değiştirilebilir.

Temel ayarlar:

* POP_SIZE: Her nesildeki birey sayısıdır. Artırılırsa daha geniş tasarım uzayı aranır ancak çözüm süresi artar.
* N_OFFSPRINGS: Her nesilde üretilecek yeni birey sayısıdır. Artırılması çeşitliliği artırabilir.
* N_EVAL: Toplam değerlendirme sayısıdır. Daha yüksek değer daha iyi Pareto dağılımı verebilir fakat analiz süresini artırır.
* SEED: Rastgelelik kontrolüdür. Aynı seed değeri aynı koşullarda benzer sonuçlar üretir. Farklı seed değerleri farklı Pareto dağılımları verebilir.
* CROSSOVER_PROB: Çaprazlama olasılığıdır. Genellikle 0.80-0.95 aralığında tutulabilir.
* CROSSOVER_ETA: Çaprazlama dağılım parametresidir. Daha yüksek değer ebeveynlere yakın bireyler üretir.
* MUTATION_ETA: Mutasyon dağılım parametresidir. Daha düşük değer daha geniş değişimlere izin verir; daha yüksek değer daha küçük değişimler üretir.

Daha hızlı deneme için:
POP_SIZE ve N_EVAL düşük tutulabilir.

Daha detaylı ve güvenilir Pareto cephesi için:
POP_SIZE ve N_EVAL artırılabilir.

ÖRNEK ANALİZ AYARLARI
POP_SIZE = 10
N_OFFSPRINGS = 5
N_EVAL = 50
SEED = 7
CROSSOVER_PROB = 0.90
CROSSOVER_ETA = 20
MUTATION_ETA = 25

NOT
Bu kodlar akademik çalışma ve bitirme projesi kapsamında hazırlanmıştır. Sonuçlar kullanılan mesh yoğunluğu, sınır koşulları, modal çözüm ayarları ve seçilen tasarım değişken aralıklarına bağlıdır.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

EN

PROJECT TITLE
Multi-Objective Structural Optimization of Stiffened Plates in Contact with Fluid

PREPARED BY
Emir Han TURAN

SUPERVISOR
Dr. Mustafa Erden YILDIZDAĞ

INSTITUTION
Istanbul Technical University
Faculty of Naval Architecture and Ocean Engineering
Department of Naval Architecture and Marine Engineering

PROJECT DESCRIPTION
This project investigates the wet modal behavior and structural mass of a stiffened plate partially in contact with fluid. The structural plate and stiffeners are modeled using SHELL181 elements, while the fluid domain is modeled using FLUID30 acoustic elements. Finite element analyses are performed in ANSYS MAPDL / PyMAPDL, and the multi-objective optimization process is carried out in Python using the pymoo library and the NSGA-II algorithm.

FILE CONTENTS

1. Solver.py
   This is the main parametric finite element solver file. It includes geometry generation, material definitions, element types, mesh settings, boundary conditions, FSI definition and wet modal analysis steps. The evaluate_design(...) function returns the selected wet natural frequency, structural volume and structural mass for a given design.

2. Test.py
   This file is used to check whether the solver works correctly. It runs a single reference analysis and prints the selected frequency, volume, mass, number of stiffeners and modal information. It is recommended to run this file before starting the optimization.

3. Optimization.py
   This is the NSGA-II optimization file. It optimizes the plate thickness, stiffener thickness, stiffener height and stiffener spacing. The objectives are to maximize the selected wet natural frequency and minimize the structural mass. At the end of the optimization, Pareto results, all evaluations and the selected optimum design are saved as CSV files.

4. MAPDL macro file
   This file is optional. It can be used in ANSYS MAPDL to visualize the model, check mesh quality and display mode shapes. The macro is mainly intended for report figures, fluid/structure mesh visualization and deformed mode shape outputs rather than for running the optimization itself.

REQUIRED LIBRARIES
The following Python libraries are required:

* numpy
* pandas
* matplotlib
* ansys-mapdl-core
* pymoo

Example installation command:
pip install numpy pandas matplotlib ansys-mapdl-core pymoo

ANSYS REQUIREMENT
ANSYS MAPDL must be installed on the computer. ANSYS 2023 R1 was used in this study. PyMAPDL is used to launch and control MAPDL analyses directly from Python.

USAGE ORDER

1. Place all Python files in the same folder.
2. Make sure that Solver.py can be imported correctly.
3. Run Test.py first to perform a single reference analysis.
4. If the test is successful, run Optimization.py to start the NSGA-II optimization.
5. After the analysis is completed, check the generated CSV result files and plot outputs in the working directory.
6. If the MAPDL macro file is used, it can be executed inside ANSYS MAPDL to obtain mesh and mode shape figures.

CHANGING NSGA-II SETTINGS
To obtain different optimization results, the NSGA-II parameters in the USER SETTINGS section of Optimization.py can be modified.

Main settings:

* POP_SIZE: Number of individuals in each generation. Increasing it explores a wider design space but increases computational time.
* N_OFFSPRINGS: Number of new individuals generated in each generation. Increasing it may improve population diversity.
* N_EVAL: Total number of function evaluations. Higher values may produce a better Pareto distribution but require longer analysis time.
* SEED: Random seed value. The same seed gives similar results under the same conditions. Different seed values may lead to different Pareto distributions.
* CROSSOVER_PROB: Crossover probability. It is commonly kept between 0.80 and 0.95.
* CROSSOVER_ETA: Crossover distribution parameter. Higher values generate offspring closer to the parent designs.
* MUTATION_ETA: Mutation distribution parameter. Lower values allow larger design changes, while higher values produce smaller changes.

For a faster trial:
Use lower POP_SIZE and N_EVAL values.

For a more detailed and reliable Pareto front:
Increase POP_SIZE and N_EVAL values.

EXAMPLE OPTIMIZATION SETTINGS
POP_SIZE = 10
N_OFFSPRINGS = 5
N_EVAL = 50
SEED = 7
CROSSOVER_PROB = 0.90
CROSSOVER_ETA = 20
MUTATION_ETA = 25

NOTE
These codes were prepared for academic use within the scope of a graduation project. The results depend on mesh density, boundary conditions, modal solution settings and the selected design variable ranges.
