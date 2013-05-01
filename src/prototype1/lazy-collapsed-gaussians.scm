(declare (usual-integrations +))
(load "randutil")

;; setting up generic extension for +

(define default+ +)
(define + (make-generic-procedure (make-procedure-arity 1 #f) 'generic+))
(set-generic-procedure-default-generator!
  + (lambda (proc tags) (force-sample-and-apply default+)))

;; sample generic datatype

(define the-unsampled-value 'the-unsampled-value)

(define sample:val (make-generic-procedure (make-procedure-arity 1) 'sample:val))
(define sample:set-val! (make-generic-procedure (make-procedure-arity 2) 'sample:set-val!))
(define sample:sampler (make-generic-procedure (make-procedure-arity 1) 'sample:sampler))
(define sample:params (make-generic-procedure (make-procedure-arity 1) 'sample:params))
(define sample? (make-generic-procedure (make-procedure-arity 1) 'sample?))

(set-generic-procedure-default-generator!
  sample? (lambda (proc tags) (lambda args #f)))

(define (sample:unsampled? s)
  (and (sample? s)
       (eq? (sample:val s) the-unsampled-value)))

(define (sample:sampled? s)
  (and (sample? s)
       (not (eq? (sample:val s) the-unsampled-value))))

(define (force-sample s)
  (if (not (sample? s))
    s
    (begin
      (if (sample:unsampled? s)
        (sample:set-val! s ((sample:sampler s) (sample:params s))))
      (sample:val s))))

(define ((force-sample-and-apply proc) . args)
  (apply proc (map force-sample args)))

(define-syntax make-sample-datatype
  (sc-macro-transformer
    (lambda (form uenv)
      (let* ((name (symbol-append (cadr form) '-sample))
             (param-maker (close-syntax (caddr form) uenv))
             (sampler (close-syntax (cadddr form) uenv))

             (record-name (symbol-append '< name '>))
             (new (symbol-append name ':new))
             (tester (symbol-append name '?))
             (param-getter (symbol-append name ':params))
             (val-getter (symbol-append name ':val))
             (val-setter (symbol-append name ':set-val!))
             (dispatch-tag-tester (symbol-append name '-dispatch-tag?)))

           `(begin
              (define-record-type ,record-name
                               (,(symbol-append '% new) params val)
                               ,tester
                               (params ,param-getter)
                               (val ,val-getter ,val-setter))

              (define (,new . args)
                (,(symbol-append '% new) (apply ,param-maker args) the-unsampled-value))

              (add-generic-procedure-generator
                sample:val (lambda (proc tags)
                             (and (eq? (length tags) 1)
                                  (,dispatch-tag-tester (car tags))
                                  ,val-getter)))

              (add-generic-procedure-generator
                sample:set-val! (lambda (proc tags)
                                  (and (eq? (length tags) 1)
                                       (,dispatch-tag-tester (car tags))
                                       ,val-setter)))

              (add-generic-procedure-generator
                sample:sampler (lambda (proc tags)
                                 (and (eq? (length tags) 1)
                                      (,dispatch-tag-tester (car tags))
                                      ,sampler)))

              (add-generic-procedure-generator
                sample? (lambda (proc tags)
                          (and (eq? (length tags) 1)
                               (,dispatch-tag-tester (car tags))
                               (lambda (s) #t))))

              (define (,dispatch-tag-tester tag)
                (eq? tag (record-type-dispatch-tag ,record-name)))
           )))))

;; gaussian

(make-sample-datatype gaussian gaussian:make-params gaussian:rvs)

;; TODO needs to handle samples (collapsed and uncollapsed, gaussian and
;; others), numbers, others are error
(define (gaussian:+ . args)
  (define (sum-means gaussians)
    (apply default+ (map gaussian:mean (map gaussian-sample:params gaussians))))
  (define (sum-vars gaussians)
    (apply default+ (map gaussian:var (map gaussian-sample:params gaussians))))
  ((partition (lambda (s) (and (gaussian-sample? s) (sample:unsampled? s))) args)
   (lambda (gaussians nongaussians)
     (let ((newmean (default+ (sum-means gaussians) (apply default+ (map force-sample nongaussians))))
           (newvar (sum-vars gaussians)))
     (gaussian-sample:new newmean newvar)))))


(add-generic-procedure-generator
  + (lambda (proc tags)
      (and (any gaussian-sample-dispatch-tag? tags)
           gaussian:+)))

